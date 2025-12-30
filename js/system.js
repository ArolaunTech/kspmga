import * as THREE from 'three';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { Orbit } from './orbit.js';
import { Vector3 } from './vector.js';
import { colorArrToHex } from './color.js';

export class System {
	bodies = [];
	bodymap = new Map();

	constructor(JSONPath, f) {
		this.ready = false;

		fetch(JSONPath)
			.then(response => response.json())
			.then(json => {
				this.bodies = json.bodies;

				let numbodies = this.bodies.length;
				for (let i = 0; i < numbodies; i++) {
					this.bodymap.set(this.bodies[i].name, i);
				}

				for (let i = 0; i < numbodies; i++) {
					if (this.bodies[i].root) continue;

					let parent = this.bodies[i].orbit.parent;
					let mu = this.bodies[this.bodymap.get(parent)].gravparameter

					let neworbit = new Orbit(
						mu, 
						this.bodies[i].orbit.sma,
						this.bodies[i].orbit.eccentricity,
						this.bodies[i].orbit.inclination,
						this.bodies[i].orbit.lan,
						this.bodies[i].orbit.argp,
						this.bodies[i].orbit.epoch,
						this.bodies[i].orbit.meananomalyatepoch,
					);

					this.bodies[i].parent = parent;
					this.bodies[i].orbit = neworbit;
				}

				f(this);

				this.ready = true;
			});
	}

	get rootmu() {
		let numbodies = this.bodies.length;
		for (let i = 0; i < numbodies; i++) {
			if (this.bodies[i].root) {
				return this.bodies[i].gravparameter;
			}
		}
		return 0;
	}

	getDState(name, time) {
		// state - state of parent
		let bodyindex = this.bodymap.get(name);
		let body = this.bodies[bodyindex];

		if (body.root) {
			return {r: new Vector3(0, 0, 0), v: new Vector3(0, 0, 0)};
		}

		return body.orbit.getState(time);
	}

	getState(name, time) {
		let bodyindex = this.bodymap.get(name);
		let body = this.bodies[bodyindex];

		if (body.root) {
			return {r: new Vector3(0, 0, 0), v: new Vector3(0, 0, 0)};
		}

		//Parent state + currrent state
		let parent = body.parent;
		let parentstate = this.getState(parent, time);
		let currstate = body.orbit.getState(time);

		return {r: Vector3.add(parentstate.r, currstate.r), v: Vector3.add(parentstate.v, currstate.v)};
	}

	fillScene(scene, scale, xres, yres) {
		let numbodies = this.bodies.length;
		let groups = [];

		const loader = new THREE.TextureLoader();
		const texture = loader.load('../sprites/circle.png');

		for (let i = 0; i < numbodies; i++) {
			groups.push(new THREE.Group());
		}

		for (let i = 0; i < numbodies; i++) {
			let colorhex = colorArrToHex(this.bodies[i].color);

			// Body sphere
			const geometry = new THREE.SphereGeometry(this.bodies[i].radius * scale, 32, 16);
			const bodymaterial = new THREE.MeshBasicMaterial({color: colorhex});
			const body = new THREE.Mesh(geometry, bodymaterial);

			groups[i].add(body);

			if (!this.bodies[i].root) {
				// SOI
				const soimaterial = new THREE.MeshBasicMaterial({color: colorhex, transparent: true, opacity: 0.5});
				const soisphere = new THREE.SphereGeometry(this.bodies[i].soi * scale, 32, 16);
				const soi = new THREE.Mesh(soisphere, soimaterial);

				groups[i].add(soi);

				// Orbit
				let parent = this.bodies[i].parent;
				let parentindex = this.bodymap.get(parent);

				let segments = 100;
				let points = [];
				for (let t = 0; t <= segments; t++) {
					let E = 2 * Math.PI * t / segments;

					let pos = this.bodies[i].orbit.getPosFromE(E)

					points.push(new THREE.Vector3(pos.x * scale, pos.z * scale, -pos.y * scale));
				}

				const orbit = new LineGeometry();
				orbit.setFromPoints(points);
				const orbitmaterial = new LineMaterial({color: colorhex, linewidth: 2});
    			orbitmaterial.resolution.set(xres, yres);
    			const orbitline = new Line2(orbit, orbitmaterial);

    			groups[parentindex].add(orbitline);

    			//Body sprite
    			if (this.bodies[parentindex].root) {
    				// Planets only
					const spritematerial = new THREE.SpriteMaterial({
						color: colorhex, 
						sizeAttenuation: false, 
						map: texture,
						depthTest: false
					});

					const sprite = new THREE.Sprite(spritematerial);
					sprite.scale.set(0.03, 0.03, 1);

					groups[i].add(sprite);
				}
			}
		}

		for (let i = 0; i < numbodies; i++) {
			scene.add(groups[i]);
		}

		return groups;
	}

	updateScene(groups, scale, time) {
		let numbodies = this.bodies.length;

		for (let i = 0; i < numbodies; i++) {
			let position = this.getState(this.bodies[i].name, time).r;

			groups[i].position.set(position.x * scale, position.z * scale, -position.y * scale);
		}
	}
};