import * as THREE from 'three';
import { Vector3 } from './vector.js';
import { Orbit, propagate } from './orbit.js';
import { lambert } from './lambert.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { colorArrToHex } from './color.js';
import { disposeHierarchy } from './dispose.js';
import { secsToKerbalTimeString } from './kerbaltime.js';
import { calcEjectiondV, calcPeriapsis, violation } from './trajectorycalcs.js';

export class Renderer {
	constructor(options) {
		this.scene = new THREE.Scene();
		this.camera = new THREE.PerspectiveCamera(
			options.fov, 
			options.aspect, 
			options.mindist, 
			options.maxdist
		);

		this.canvas = options.htmlCanvas;
		this.renderer = new THREE.WebGLRenderer({canvas: this.canvas, antialias: true});
		this.renderer.setSize(this.canvas.width, this.canvas.height, false);
		this.renderer.setAnimationLoop(this.createAnimationCallback(this));

		this.controls = new OrbitControls(this.camera, this.canvas);
		this.controls.minDistance = 0.5;
		this.controls.maxDistance = 500;
		this.camera.position.z = 5;
		this.controls.update();

		this.scale = options.scale;

		this.groups = [];
		this.trajectorygroup = new THREE.Group();
		this.scene.add(this.trajectorygroup);

		this.loader = new THREE.TextureLoader();

		this.podsprite = this.createPodSprite();
	}

	createPodSprite() {
		const texture = this.loader.load("https://arolauntech.github.io/kspmga/sprites/pod.png");

		const spritematerial = new THREE.SpriteMaterial({
			sizeAttenuation: false, 
			map: texture,
			depthTest: false
		});

		const sprite = new THREE.Sprite(spritematerial);
		sprite.scale.set(0.05, 0.05, 1);
		sprite.visible = false;

		this.scene.add(sprite);

		return sprite;
	}

	fillSceneWithSystem(system) {
		let numbodies = system.bodies.length;
		this.groups = [];

		const texture = this.loader.load('https://arolauntech.github.io/kspmga/sprites/circle.png');

		for (let i = 0; i < numbodies; i++) {
			this.groups.push(new THREE.Group());
		}

		for (let i = 0; i < numbodies; i++) {
			let colorhex = colorArrToHex(system.bodies[i].color);

			// Body sphere
			const geometry = new THREE.SphereGeometry(system.bodies[i].radius * this.scale, 32, 16);
			const bodymaterial = new THREE.MeshBasicMaterial({color: colorhex});
			const body = new THREE.Mesh(geometry, bodymaterial);

			this.groups[i].add(body);

			if (!system.bodies[i].root) {
				// SOI
				const soimaterial = new THREE.MeshBasicMaterial({color: colorhex, transparent: true, opacity: 0.5});
				const soisphere = new THREE.SphereGeometry(system.bodies[i].soi * this.scale, 32, 16);
				const soi = new THREE.Mesh(soisphere, soimaterial);

				this.groups[i].add(soi);

				// Orbit
				let parent = system.bodies[i].parent;
				let parentindex = system.bodymap.get(parent);

				let segments = 100;
				let points = [];
				for (let t = 0; t <= segments; t++) {
					let E = 2 * Math.PI * t / segments;

					let pos = system.bodies[i].orbit.getPosFromE(E);

					points.push(new THREE.Vector3(pos.x * this.scale, pos.z * this.scale, -pos.y * this.scale));
				}

				const orbit = new LineGeometry();
				orbit.setFromPoints(points);
				const orbitmaterial = new LineMaterial({color: colorhex, linewidth: 2});
				orbitmaterial.resolution.set(this.canvas.width, this.canvas.height);
				const orbitline = new Line2(orbit, orbitmaterial);

				this.groups[parentindex].add(orbitline);

				//Body sprite
				if (system.bodies[parentindex].root) {
					// Planets only
					const spritematerial = new THREE.SpriteMaterial({
						color: colorhex, 
						sizeAttenuation: false, 
						map: texture,
						depthTest: false
					});

					const sprite = new THREE.Sprite(spritematerial);
					sprite.scale.set(0.03, 0.03, 1);

					this.groups[i].add(sprite);
				}
			}
		}

		for (let i = 0; i < numbodies; i++) {
			this.scene.add(this.groups[i]);
		}
	}

	updateSceneWithSystem(system, time) {
		let numbodies = system.bodies.length;

		for (let i = 0; i < numbodies; i++) {
			let position = system.getState(system.bodies[i].name, time).r;

			this.groups[i].position.set(position.x * this.scale, position.z * this.scale, -position.y * this.scale);
		}
	}

	displayTrajectory(trajectory, system) {
		console.log(trajectory);
		console.log(JSON.stringify(trajectory));

		disposeHierarchy(this.trajectorygroup);
		this.trajectorygroup.clear();

		this.podsprite.visible = true;

		let sequence = trajectory[2][0];
		let parent = system.bodies[system.bodymap.get(sequence[0])].parent;
		let mu = system.bodies[system.bodymap.get(parent)].gravparameter;

		let trajectorydetailstring = `Optimized MGA trajectory details:\n\n${sequence[0]} departure:\n`;
		trajectorydetailstring += ` - Ejection Δv: ${calcEjectiondV(sequence[0], Math.hypot(trajectory[0][1], trajectory[0][2], trajectory[0][3]), trajectory[2][1], system).toFixed(3)} m/s\n`;
		trajectorydetailstring += ` - Time: ${secsToKerbalTimeString(trajectory[0][0])}\n`;

		for (let i = 0; i < sequence.length - 1; i++) {
			let t1 = trajectory[0][5 * i];
			let t2 = trajectory[0][5 * i + 5];

			let state1 = system.getDState(sequence[i], t1);
			let state2 = system.getDState(sequence[i + 1], t2);

			let outgoing = new Vector3(
				trajectory[0][5 * i + 1],
				trajectory[0][5 * i + 2],
				trajectory[0][5 * i + 3]
			);

			let vout = Vector3.add(outgoing, state1.v);

			let todsm = trajectory[0][5 * i + 4] * (t2 - t1);
			let tdsm = t1 + todsm;

			let sma = 1 / (2 / state1.r - (vout.norm * vout.norm) / mu);
			let period = Infinity;
			if (sma > 0) {
				period = 2 * Math.PI * Math.sqrt(sma * sma * sma / mu);
			}

			if (todsm > period) todsm = period;

			const segments = 100;
			let points = [];
			for (let j = 0; j <= segments; j++) {
				let time = todsm * j / segments;

				let interstate = propagate(state1.r, vout, time, mu);

				points.push(new THREE.Vector3(interstate.r.x * this.scale, interstate.r.z * this.scale, -interstate.r.y * this.scale));
			}

			const outorbitdisplay = new LineGeometry();
			outorbitdisplay.setFromPoints(points);
			const outorbitmaterial = new LineMaterial({color: 0xffffff, linewidth: 2});
			outorbitmaterial.resolution.set(this.canvas.width, this.canvas.height);
			const outorbitline = new Line2(outorbitdisplay, outorbitmaterial);

			this.trajectorygroup.add(outorbitline);

			let statedsm = propagate(state1.r, vout, trajectory[0][5 * i + 4] * (t2 - t1), mu);
			let sols = lambert(statedsm.r, state2.r, (1 - trajectory[0][5 * i + 4]) * (t2 - t1), mu, 5, false);

			let mindvdsm = Infinity;
			let bestvafter;
			let bestvin;

			for (let j = 0; j < sols.length; j++) {
				let dvdsm = Vector3.sub(sols[j][0], statedsm.v).norm;

				if (dvdsm < mindvdsm) {
					mindvdsm = dvdsm;
					bestvafter = sols[j][0];
					bestvin = Vector3.sub(sols[j][1], state2.v);
				}
			}

			let fromdsm = (1 - trajectory[0][5 * i + 4]) * (t2 - t1);

			sma = 1 / (2 / statedsm.r - (bestvafter.norm * bestvafter.norm) / mu);
			period = Infinity;
			if (sma > 0) {
				period = 2 * Math.PI * Math.sqrt(sma * sma * sma / mu);
			}

			if (fromdsm > period) fromdsm = period;

			points = [];
			for (let j = 0; j <= segments; j++) {
				let time = fromdsm * j / segments;

				let interstate = propagate(statedsm.r, bestvafter, time, mu);

				points.push(new THREE.Vector3(interstate.r.x * this.scale, interstate.r.z * this.scale, -interstate.r.y * this.scale));
			}

			const inorbitdisplay = new LineGeometry();
			inorbitdisplay.setFromPoints(points);
			const inorbitmaterial = new LineMaterial({color: 0xff0000, linewidth: 2});
			inorbitmaterial.resolution.set(this.canvas.width, this.canvas.height);
			const inorbitline = new Line2(inorbitdisplay, inorbitmaterial);

			this.trajectorygroup.add(inorbitline);

			let progradedirdsm = Vector3.normalize(statedsm.v);
			let normaldirdsm = Vector3.normalize(Vector3.cross(statedsm.r, statedsm.v));
			let radialdirdsm = Vector3.normalize(Vector3.cross(progradedirdsm, normaldirdsm));

			trajectorydetailstring += `\n${sequence[i]}-${sequence[i + 1]} DSM:\n`;
			trajectorydetailstring += ` - Δv: ${mindvdsm.toFixed(3)} m/s\n`;
			trajectorydetailstring += ` - Time: ${secsToKerbalTimeString(tdsm)}\n`;
			trajectorydetailstring += ` - Prograde: ${Vector3.dot(Vector3.sub(bestvafter, statedsm.v), progradedirdsm).toFixed(3)} m/s\n`;
			trajectorydetailstring += ` - Normal: ${Vector3.dot(Vector3.sub(bestvafter, statedsm.v), normaldirdsm).toFixed(3)} m/s\n`;
			trajectorydetailstring += ` - Radial out: ${Vector3.dot(Vector3.sub(bestvafter, statedsm.v), radialdirdsm).toFixed(3)} m/s\n`;

			if ((i === sequence.length - 2) && (trajectory[2][10])) {
				trajectorydetailstring += `\n${sequence[i + 1]} capture:\n`;
			} else {
				trajectorydetailstring += `\n${sequence[i + 1]} flyby:\n`;
			}

			trajectorydetailstring += ` - Time: ${secsToKerbalTimeString(t2)}\n`;
			trajectorydetailstring += ` - Incoming v-inf: ${bestvin.norm.toFixed(3)} m/s\n`;
			if (i === sequence.length - 2) {
				trajectorydetailstring += ` - Periapsis: ${(trajectory[2][2] / 1000).toFixed(3)} km\n`;
				trajectorydetailstring += ` - Required insertion: ${calcEjectiondV(sequence[sequence.length - 1], bestvin.norm, trajectory[2][2], system).toFixed(3)} m/s\n`;
			} else {
				let nextvout = new Vector3(trajectory[0][5 * i + 6], trajectory[0][5 * i + 7], trajectory[0][5 * i + 8]);

				trajectorydetailstring += ` - Periapsis: ${(calcPeriapsis(sequence[i + 1], bestvin, nextvout, system) / 1000).toFixed(3)} km\n`;
				trajectorydetailstring += ` - Flyby burn at exit: ${(violation(sequence[i + 1], bestvin, nextvout, system)).toFixed(3)} m/s\n`;
			}
		}

		return trajectorydetailstring;
	}

	updatePodTrajectory(trajectory, system, time) {
		let sequence = trajectory[2][0];
		let parent = system.bodies[system.bodymap.get(sequence[0])].parent;
		let mu = system.bodies[system.bodymap.get(parent)].gravparameter;

		for (let i = 0; i < sequence.length - 1; i++) {
			let t1 = trajectory[0][5 * i];
			let t2 = trajectory[0][5 * i + 5];

			if ((time > t2) && (i < sequence.length - 2)) continue;

			let state1 = system.getDState(sequence[i], t1);
			let state2 = system.getDState(sequence[i + 1], t2);

			let outgoing = new Vector3(
				trajectory[0][5 * i + 1],
				trajectory[0][5 * i + 2],
				trajectory[0][5 * i + 3]
			);

			let vout = Vector3.add(outgoing, state1.v);

			let todsm = trajectory[0][5 * i + 4] * (t2 - t1);
			let tdsm = t1 + todsm;

			if (time < tdsm) {
				let podstate = propagate(state1.r, vout, (time - t1), mu);

				this.podsprite.position.set(podstate.r.x * this.scale, podstate.r.z * this.scale, -podstate.r.y * this.scale);
				return;
			}

			let statedsm = propagate(state1.r, vout, todsm, mu);

			let sols = lambert(statedsm.r, state2.r, (1 - trajectory[0][5 * i + 4]) * (t2 - t1), mu, 5, false);

			let mindvdsm = Infinity;
			let bestvafter;
			let bestvin;

			for (let j = 0; j < sols.length; j++) {
				let dvdsm = Vector3.sub(sols[j][0], statedsm.v).norm;

				if (dvdsm < mindvdsm) {
					mindvdsm = dvdsm;
					bestvafter = sols[j][0];
					bestvin = Vector3.sub(sols[j][1], state2.v);
				}
			}

			let podstate = propagate(statedsm.r, bestvafter, time - tdsm, mu);

			this.podsprite.position.set(podstate.r.x * this.scale, podstate.r.z * this.scale, -podstate.r.y * this.scale);
			
			return;
		}
	}

	createAnimationCallback(renderer) {
		return function() {
			renderer.controls.update();
			renderer.renderer.setSize(renderer.canvas.width, renderer.canvas.height, false);
			renderer.renderer.render(renderer.scene, renderer.camera);
		};
	}
};