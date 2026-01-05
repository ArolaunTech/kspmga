import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { System } from './system.js';
import { fillScene, updateScene } from './scenefiller.js';
import { secsToKerbalTimeString, kerbalTimeToSecs } from './kerbaltime.js';

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

const initbody = document.getElementById("startingbody");
const finalbody = document.getElementById("endingbody");
const initalt = document.getElementById("initalt");
const finalalt = document.getElementById("finalalt");
const minvinf = document.getElementById("minvinf");
const maxvinf = document.getElementById("maxvinf");
const maxdvdsm = document.getElementById("maxdvdsm");
const maxduration = document.getElementById("maxduration");
const earliestyear = document.getElementById("earliestyear");
const earliestday = document.getElementById("earliestday");
const earliesthour = document.getElementById("earliesthour");
const latestyear = document.getElementById("latestyear");
const latestday = document.getElementById("latestday");
const latesthour = document.getElementById("latesthour");
const includeinsertion = document.getElementById("includecapture");

const startsearch = document.getElementById("startsearch");

const errormsg = document.getElementById("errormsg");

function handleNumericMinMax() {
	if (this.value === "") return;

	let numericvalue = Number(this.value);

	if ((this.min !== "") && (numericvalue < Number(this.min))) {
		this.value = this.min;
	}

	if ((this.max !== "") && (numericvalue > Number(this.max))) {
		this.value = this.max;
	}
}

function handleNumericMinMaxInt() {
	if (this.value === "") return;

	let numericvalue = Number(this.value);

	if ((this.min !== "") && (numericvalue < Number(this.min))) {
		this.value = this.min;
	}

	if ((this.max !== "") && (numericvalue > Number(this.max))) {
		this.value = this.max;
	}

	if (!Number.isInteger(numericvalue)) {
		this.value = Math.round(numericvalue);
	}
}

function setBodySelectOptions(element, sys) {
	element.textContent = "";

	for (let i = 0; i < sys.bodies.length; i++) {
		if (sys.bodies[i].root) continue; // Cannot currently do flyby sequences with arbitrary solar orbits

		const option = document.createElement("option");

		option.innerText = sys.bodies[i].name;

		element.appendChild(option);
	}
}

[initalt, finalalt, minvinf, maxvinf, maxdvdsm, maxduration].forEach((element) => {
	element.onchange = handleNumericMinMax;
});

[earliestyear, earliesthour, earliestday, latestyear, latesthour, latestday].forEach((element) => {
	element.onchange = handleNumericMinMaxInt;
});

timeslider.oninput = function() {
	timedisplay.innerText = secsToKerbalTimeString(this.value);
	time = this.value;
	needsupdate = true;
}

function updateCanvasSize() {
	if (document.body.clientWidth > 1000) {
		canvas.width = 640;
	} else {
		canvas.width = document.body.clientWidth * 0.64;
	}
	canvas.height = canvas.width / aspectratio;
}

window.onresize = updateCanvasSize;

startsearch.onclick = function() {
	errormsg.innerText = "";

	let valid = sequenceregex.test(sequence.value);

	if (!valid) {
		errormsg.innerText = "Error: invalid sequence provided!";
		return;
	}

	if (initalt.value.length === 0) {
		errormsg.innerText = "Error: must provide an initial altitude!";
		return;
	}

	if (finalalt.value.length === 0) {
		errormsg.innerText = "Error: must provide a final altitude!";
		return;
	}

	let initaltnum = 1000 * Number(initalt.value);
	let finalaltnum = 1000 * Number(finalalt.value);

	if (minvinf.value.length === 0) {
		errormsg.innerText = "Error: must provide a minimum v-inf!";
		return;
	}

	if (maxvinf.value.length === 0) {
		errormsg.innerText = "Error: must provide a maximum v-inf!";
		return;
	}

	let minvinfnum = Number(minvinf.value);
	let maxvinfnum = Number(maxvinf.value);

	if (maxvinfnum < minvinfnum) {
		errormsg.innerText = "Error: maximum v-inf must be larger than minimum v-inf!";
		return;
	}

	if (maxdvdsm.value.length === 0) {
		errormsg.innerText = "Error: must provide a maximum Δv per DSM!";
		return;
	}

	let maxdvdsmnum = Number(maxdvdsm.value);

	let maxdurationnum = Number(maxduration.value);

	if (maxduration.value.length === 0) {
		maxdurationnum = Infinity;
	}

	if ((earliestyear.value.length === 0) || (earliestday.value.length === 0) || (earliesthour.value.length === 0)) {
		errormsg.innerText = "Error: must provide an initial date!";
		return;
	}

	let earliesttime = kerbalTimeToSecs([Number(earliestyear.value), Number(earliestday.value), Number(earliesthour.value), 0, 0]);
	let latesttime = kerbalTimeToSecs([Number(latestyear.value), Number(latestday.value), Number(latesthour.value), 0, 0]);

	if ((latestyear.value.length === 0) || (latestday.value.length === 0) || (latesthour.value.length === 0)) {
		latesttime = Infinity;
	}

	console.log("h", mgafinder);

	mgafinder.postMessage({
		init: false,
		params: [
			sequence.value.split("-"),
			initaltnum,
			finalaltnum,
			minvinfnum, 
			maxvinfnum,
			maxdvdsmnum,
			maxdurationnum,
			maxrevsnum,
			earliesttime,
			latesttime,
			includecapture.checked,
			includeflyby.checked
		]
	});
}

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
	groups = fillScene(sys, scene, SCALE, canvas.width, canvas.height);

	setBodySelectOptions(initbody, sys);
	setBodySelectOptions(finalbody, sys);

	mgafinder = new Worker(new URL("mga.js", import.meta.url), {type: "module"});

	mgafinder.postMessage({init: true, system: sys});
});

camera.position.z = 5;
controls.update();

function animate() {
	controls.update();
	renderer.setSize(canvas.width, canvas.height, false);
	renderer.render(scene, camera);

	if (system.ready && needsupdate) {
		updateScene(system, groups, SCALE, time);
		needsupdate = false;
	}
}