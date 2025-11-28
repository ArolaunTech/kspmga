import * as THREE from 'three';

let system = new System("../data/systems/stock.json");

console.log(system);

const canvas = document.getElementById("system");

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, canvas.width / canvas.height, 0.1, 1000);

const renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
renderer.setSize(canvas.width, canvas.height);
renderer.setAnimationLoop(animate);

const geometry = new THREE.SphereGeometry(1, 32, 16);
const material = new THREE.MeshStandardMaterial({
	color: 0x00ff00,
	emissive: 0xffff80
});
const sphere = new THREE.Mesh(geometry, material);
scene.add(sphere);

camera.position.z = 5;

function animate() {
	sphere.rotation.x += 0.01;
	sphere.rotation.y += 0.01;

	renderer.render(scene, camera);
}