import * as THREE from 'three';
import { Vector3 } from './vector.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { System } from './system.js';
import { fillScene, updateScene } from './scenefiller.js';
import { secsToKerbalTimeString, kerbalTimeToSecs } from './kerbaltime.js';
import { handleNumericMinMax, handleNumericMinMaxInt, setBodySelectOptions } from './input.js';

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

let disabled = false;

const addbody = document.getElementById("addbody");
const subbody = document.getElementById("subbody");

let numbodies = 2;

let intermediatebodies = [];
let flybydvs = [];
let maxrevs = [document.getElementById("maxrevs")];
let inputstack = [];

const startsearch = document.getElementById("startsearch");

const errormsg = document.getElementById("errormsg");

let trajectorygroup = new THREE.Group();
function displayTrajectory(trajectory) {
	console.log(trajectory);
	console.log(JSON.stringify(trajectory));

	

	let sequence = trajectory[2][0];
	for (let i = 0; i < sequence.length - 1; i++) {
		let t1 = trajectory[0][5 * i];
		let state1 = system.getDState(sequence[i], t1);
		
		let outgoing = new Vector3(
			trajectory[0][5 * i + 1],
			trajectory[0][5 * i + 2],
			trajectory[0][5 * i + 3]
		);

		let vout = Vector3.add(outgoing, state1.v);
	}
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

timeslider.oninput = function() {
	timedisplay.innerText = secsToKerbalTimeString(this.value);
	time = this.value;
	needsupdate = true;
}

startsearch.onclick = function() {
	errormsg.innerText = "";

	if (disabled) return;

	let sequence = [initbody.value];

	for (let i = 0; i < intermediatebodies.length; i++) {
		sequence.push(intermediatebodies[i].value);
	}

	sequence.push(finalbody.value);

	console.log(sequence);

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

	let maxdurationnum = Number(maxduration.value) * 21600;

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

	let flybydvnums = [];
	for (let i = 0; i < flybydvs.length; i++) {
		if (flybydvs[i].value.length === 0) {
			errormsg.innerText = "Error: must provide a maximum flyby Δv!";
			return;
		}

		flybydvnums.push(Number(flybydvs[i].value));
	}

	let maxrevnums = [];
	for (let i = 0; i < maxrevs.length; i++) {
		if (maxrevs[i].value.length === 0) {
			errormsg.innerText = "Error: must provide a maximum number of revolutions for all transfers!";
			return;
		}

		maxrevnums.push(Number(maxrevs[i].value));
	}

	console.log("h", mgafinder);

	mgafinder.postMessage({
		init: false,
		params: [
			sequence,
			initaltnum,
			finalaltnum,
			minvinfnum, 
			maxvinfnum,
			maxdvdsmnum,
			maxdurationnum,
			maxrevnums,
			earliesttime,
			latesttime,
			includeinsertion.checked,
			flybydvnums
		]
	});

	disabled = true;
}

addbody.onclick = function() {
	const flybysettings = document.createElement("div");
	const descriptor = document.createElement("p");

	flybysettings.classList.add("startsearchbox");

	descriptor.innerText = "Flyby";

	const bodygroup = document.createElement("div");
	const bodylabel = document.createElement("label");
	const bodyinput = document.createElement("select");

	bodygroup.classList.add("controlgroup");

	bodylabel.classList.add("controllabel");
	bodylabel.innerText = "Flyby body";

	bodygroup.appendChild(bodylabel);

	setBodySelectOptions(bodyinput, system);

	bodygroup.appendChild(bodyinput);

	const flybydvgroup = document.createElement("div");
	const flybydvlabel = document.createElement("label");
	const flybydvinput = document.createElement("input");
	const flybydvunits = document.createElement("label");

	flybydvgroup.classList.add("controlgroup");

	flybydvlabel.classList.add("controllabel");
	flybydvlabel.innerText = "Max flyby Δv";

	flybydvinput.type = "number";
	flybydvinput.min = "0";
	flybydvinput.step = "any";
	flybydvinput.required = true;
	flybydvinput.value = "100";
	flybydvinput.onchange = handleNumericMinMax;

	flybydvunits.innerText = "m/s";

	const transferdescriptor = document.createElement("p");

	transferdescriptor.innerText = "Transfer";

	const transfersettings = document.createElement("div");
	const maxrevslabel = document.createElement("label");
	const maxrevsinput = document.createElement("input");

	transfersettings.classList.add("controlgroup");

	maxrevslabel.classList.add("controllabel");
	maxrevslabel.innerText = "Max. revolutions";

	maxrevsinput.type = "number";
	maxrevsinput.min = "0";
	maxrevsinput.step = "1";
	maxrevsinput.max = "20";
	maxrevsinput.required = true;
	maxrevsinput.value = "1";
	maxrevsinput.onchange = handleNumericMinMaxInt;

	flybydvgroup.appendChild(flybydvlabel);
	flybydvgroup.appendChild(flybydvinput);
	flybydvgroup.appendChild(flybydvunits);

	transfersettings.appendChild(maxrevslabel);
	transfersettings.appendChild(maxrevsinput);

	flybysettings.appendChild(descriptor);
	flybysettings.appendChild(bodygroup);
	flybysettings.appendChild(flybydvgroup);

	const break1 = document.createElement("br");
	const break2 = document.createElement("br");

	inputstack.push([break1, flybysettings, break2, transferdescriptor, transfersettings]);

	document.getElementById("inputs").insertBefore(break1, document.getElementById("endsettings"));
	document.getElementById("inputs").insertBefore(flybysettings, document.getElementById("endsettings"));
	document.getElementById("inputs").insertBefore(break2, document.getElementById("endsettings"));
	document.getElementById("inputs").insertBefore(transferdescriptor, document.getElementById("endsettings"));
	document.getElementById("inputs").insertBefore(transfersettings, document.getElementById("endsettings"));

	numbodies++;

	intermediatebodies.push(bodyinput);
	flybydvs.push(flybydvinput);
	maxrevs.push(maxrevsinput);
}

subbody.onclick = function() {
	if (numbodies <= 2) return;

	let removeelems = inputstack.pop();

	for (let i = 0; i < removeelems.length; i++) {
		removeelems[i].remove();
	}

	maxrevs.pop();
	flybydvs.pop();
	intermediatebodies.pop();

	numbodies--;
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
	mgafinder.onmessage = function(e) {
		if (e.data.status === "failure") return;

		disabled = false;

		displayTrajectory(e.data.result);
	}

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