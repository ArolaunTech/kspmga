import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { System } from './system.js';

const SCALE = 1e-9;
const mindist = 0.1;
const maxdist = 1000;

const canvas = document.getElementById("system");

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, canvas.width / canvas.height, mindist, maxdist);
const controls = new OrbitControls(camera, canvas);
controls.minDistance = 0.5;
controls.maxDistance = 500;

const renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
renderer.setSize(canvas.width, canvas.height);
renderer.setAnimationLoop(animate);

let groups = [];
let system = new System("../data/systems/stock.json", (sys) => {
	groups = sys.fillScene(scene, SCALE, canvas.width, canvas.height);
});

camera.position.z = 5;
controls.update();

let time = 0;
function animate() {
	controls.update();
	renderer.render(scene, camera);

	if (system.ready) {
		system.updateScene(groups, SCALE, time);
		time += 2000;
	}
}