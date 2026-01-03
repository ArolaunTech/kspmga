import * as THREE from 'three';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { colorArrToHex } from './color.js';

export function fillScene(system, scene, scale, xres, yres) {
	let numbodies = system.bodies.length;
	let groups = [];

	const loader = new THREE.TextureLoader();
	const texture = loader.load('https://arolauntech.github.io/kspmga/sprites/circle.png');

	for (let i = 0; i < numbodies; i++) {
		groups.push(new THREE.Group());
	}

	for (let i = 0; i < numbodies; i++) {
		let colorhex = colorArrToHex(system.bodies[i].color);

		// Body sphere
		const geometry = new THREE.SphereGeometry(system.bodies[i].radius * scale, 32, 16);
		const bodymaterial = new THREE.MeshBasicMaterial({color: colorhex});
		const body = new THREE.Mesh(geometry, bodymaterial);

		groups[i].add(body);

		if (!system.bodies[i].root) {
			// SOI
			const soimaterial = new THREE.MeshBasicMaterial({color: colorhex, transparent: true, opacity: 0.5});
			const soisphere = new THREE.SphereGeometry(system.bodies[i].soi * scale, 32, 16);
			const soi = new THREE.Mesh(soisphere, soimaterial);

			groups[i].add(soi);

			// Orbit
			let parent = system.bodies[i].parent;
			let parentindex = system.bodymap.get(parent);

			let segments = 100;
			let points = [];
			for (let t = 0; t <= segments; t++) {
				let E = 2 * Math.PI * t / segments;

				let pos = system.bodies[i].orbit.getPosFromE(E)

				points.push(new THREE.Vector3(pos.x * scale, pos.z * scale, -pos.y * scale));
			}

			const orbit = new LineGeometry();
			orbit.setFromPoints(points);
			const orbitmaterial = new LineMaterial({color: colorhex, linewidth: 2});
			orbitmaterial.resolution.set(xres, yres);
			const orbitline = new Line2(orbit, orbitmaterial);

			groups[parentindex].add(orbitline);

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

				groups[i].add(sprite);
			}
		}
	}

	for (let i = 0; i < numbodies; i++) {
		scene.add(groups[i]);
	}

	return groups;
}

export function updateScene(system, groups, scale, time) {
	let numbodies = system.bodies.length;

	for (let i = 0; i < numbodies; i++) {
		let position = system.getState(system.bodies[i].name, time).r;

		groups[i].position.set(position.x * scale, position.z * scale, -position.y * scale);
	}
}