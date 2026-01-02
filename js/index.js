import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { System } from './system.js';
import { MGAFinder } from './mga.js';
import { secsToKerbalTimeString } from './kerbaltime.js';

// Consts
const aspectratio = 4 / 3;
const fov = 75;
const SCALE = 1e-9;
const mindist = 0.1;
const maxdist = 1000;

// Inputs
let time = 0;
let needsupdate = true;

const timeslider = document.getElementById("timeslider");
const timedisplay = document.getElementById("timelabel");
const canvas = document.getElementById("system");

timeslider.oninput = function() {
	timedisplay.innerText = secsToKerbalTimeString(this.value);
	time = this.value;
	needsupdate = true;
}

function updateCanvasSize() {
	if (document.body.clientWidth > 1000) {
		canvas.width = 640;
	} else if (document.body.clientWidth > 500) {
		canvas.width = document.body.clientWidth * 0.64;
	} else {
		canvas.width = document.body.clientWidth - 180;
	}
	canvas.height = canvas.width / aspectratio;
}

window.onresize = updateCanvasSize;

updateCanvasSize();

// Render
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(fov, aspectratio, mindist, maxdist);
const controls = new OrbitControls(camera, canvas);
controls.minDistance = 0.5;
controls.maxDistance = 500;

const renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
renderer.setSize(canvas.width, canvas.height, false);
renderer.setAnimationLoop(animate);

let groups = [];
let mgafinder;
let system = new System("https://arolauntech.github.io/kspmga/data/systems/stock.json", (sys) => {
	groups = sys.fillScene(scene, SCALE, canvas.width, canvas.height);
	mgafinder = new MGAFinder(sys);

	// Test: Kerbin to Eve starting at 12636864 seconds with a rel. vel of 960 m/s
	let start = performance.now();
	const runs = 1;

	for (let i = 0; i < runs; i++) {
		console.log(mgafinder.findTransfersNoDSM("Kerbin", "Eve", 12636864, 960));
	}
	console.log((performance.now() - start) / runs);
});

camera.position.z = 5;
controls.update();

function animate() {
	controls.update();
	renderer.setSize(canvas.width, canvas.height, false);
	renderer.render(scene, camera);

	if (system.ready && needsupdate) {
		system.updateScene(groups, SCALE, time);
		needsupdate = false;
	}
}