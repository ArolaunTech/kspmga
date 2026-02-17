import { System } from './system.js';
import { secsToTimeString, timeToSecs, setTimeSystem, getTimeSystem } from './kerbaltime.js';
import { handleNumericMinMax, handleNumericMinMaxInt, setBodySelectOptions } from './input.js';
import { Renderer } from './renderer.js';

import { Vector3 } from './vector.js';
import { calcEjectionDetailsInclined } from './trajectorycalcs.js';

// Consts
const aspectratio = 4 / 3;

// Inputs
let time = 0;

const canvas = document.getElementById("system");
const finalalt = document.getElementById("finalalt");
const finalbody = document.getElementById("endingbody");
const initalt = document.getElementById("initalt");
const initbody = document.getElementById("startingbody");
const maxduration = document.getElementById("maxduration");
const maxdvdsm = document.getElementById("maxdvdsm");
const maxvinf = document.getElementById("maxvinf");
const minvinf = document.getElementById("minvinf");
const timeslider = document.getElementById("timeslider");
const timedisplay = document.getElementById("timelabel");
const earliestyear = document.getElementById("earliestyear");
const earliestday = document.getElementById("earliestday");
const earliesthour = document.getElementById("earliesthour");
const latestyear = document.getElementById("latestyear");
const latestday = document.getElementById("latestday");
const latesthour = document.getElementById("latesthour");
const includeinsertion = document.getElementById("includecapture");
const addbody = document.getElementById("addbody");
const subbody = document.getElementById("subbody");
const startsearch = document.getElementById("startsearch");
const stopsearch = document.getElementById("stopsearch");
const errormsg = document.getElementById("errormsg");
const trajdetails = document.getElementById("trajdetailsp");
const searchdetails = document.getElementById("searchdetails");
const systemselect = document.getElementById("systemselect");

let disabled = false;

let numbodies = 2;
let intermediatebodies = [];
let flybydvs = [];
let maxrevs = [document.getElementById("maxrevs")];
let inputstack = [];

let displayedtrajectory = [];

function updateCanvasSize() {
	if (document.body.clientWidth > 1000) {
		canvas.width = 640;
	} else {
		canvas.width = document.body.clientWidth * 0.64;
	}
	canvas.height = canvas.width / aspectratio;
}

window.onresize = updateCanvasSize;

systemselect.onchange = function() {
	system = new System(`https://arolauntech.github.io/kspmga/data/systems/${systemselect.value}.json`, (sys, json) => {
		renderer.updateSettings({
			scale: json.SCALE,
			minscrolldist: json.mindist,
			maxscrolldist: json.maxdist,
		});

		displayedtrajectory = [];

		renderer.clearTrajectory();
		renderer.fillSceneWithSystem(sys);
		renderer.updateSceneWithSystem(sys, 0);

		setBodySelectOptions(initbody, sys);
		setBodySelectOptions(finalbody, sys);

		setTimeSystem(json.timesystem);

		timeslider.min = 0;
		timeslider.max = ((getTimeSystem() == "KERBAL") ? 9201600 : 31536000);

		timedisplay.innerText = secsToTimeString(timeslider.value);
		time = timeslider.value;
	});
}

timeslider.oninput = function() {
	timedisplay.innerText = secsToTimeString(this.value);
	time = this.value;
	renderer.updateSceneWithSystem(system, time);

	if (displayedtrajectory.length > 0) renderer.updatePodTrajectory(displayedtrajectory, system, time);
}

function startMGASearch() {
	errormsg.innerText = "";

	if (disabled) return;

	let sequence = [initbody.value];
	let parent = system.bodies[system.bodymap.get(initbody.value)].parent;
	if (!system.bodies[system.bodymap.get(parent)].root) {
		errormsg.innerText = "This calculator only supports planets currently.";
		return;
	}

	for (let i = 0; i < intermediatebodies.length; i++) {
		sequence.push(intermediatebodies[i].value);

		parent = system.bodies[system.bodymap.get(sequence[i + 1])].parent;
		if (!system.bodies[system.bodymap.get(parent)].root) {
			errormsg.innerText = "This calculator only supports planets currently.";
			return;
		}
	}

	sequence.push(finalbody.value);
	parent = system.bodies[system.bodymap.get(finalbody.value)].parent;
	if (!system.bodies[system.bodymap.get(parent)].root) {
		errormsg.innerText = "This calculator only supports planets currently.";
		return;
	}

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

	let maxdurationnum = Number(maxduration.value) * ((getTimeSystem() == "KERBAL") ? 21600 : 86400);

	if (maxduration.value.length === 0) {
		maxdurationnum = Infinity;
	}

	if ((earliestyear.value.length === 0) || (earliestday.value.length === 0) || (earliesthour.value.length === 0)) {
		errormsg.innerText = "Error: must provide an initial date!";
		return;
	}

	if ((getTimeSystem() === "REAL") && ((Number(earliestyear.value) < 1951) || (Number(latestyear.value) < 1951))) {
		errormsg.innerText = "Error: must provide a valid date!";
		return;
	}

	let earliesttime = timeToSecs([Number(earliestyear.value), Number(earliestday.value), Number(earliesthour.value), 0, 0]);
	let latesttime = timeToSecs([Number(latestyear.value), Number(latestday.value), Number(latesthour.value), 0, 0]);

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

	mgafinder = new Worker(new URL("mga.js", import.meta.url), {type: "module"});
	mgafinder.onmessage = function(e) {
		if (e.data.status !== "report") {
			disabled = false;
		} else {
			searchdetails.innerText = `Found ${e.data.foundtrajectories} trajectories - ${e.data.steps} steps in`;
			return;
		}

		if (e.data.status === "failure") return;

		displayedtrajectory = e.data.result;

		trajdetailsp.innerText = renderer.displayTrajectory(e.data.result, system);

		timeslider.min = e.data.result[0][0];
		timeslider.max = e.data.result[0][e.data.result[0].length - 1];

		timedisplay.innerText = secsToTimeString(timeslider.value);
		time = timeslider.value;
		renderer.updateSceneWithSystem(system, time);
		renderer.updatePodTrajectory(e.data.result, system, time);

		mgafinder.onmessage = undefined;
		mgafinder.terminate();
		mgafinder = undefined;
	}

	mgafinder.postMessage({
		system: system,
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

function stopMGASearch() {
	mgafinder.onmessage = undefined;
	mgafinder.terminate();
	mgafinder = undefined;

	disabled = false;
}

startsearch.onclick = startMGASearch;
stopsearch.onclick = stopMGASearch;

addbody.onclick = function() {
	const flybysettings = document.createElement("div");
	flybysettings.classList.add("startsearchbox");

	const descriptor = document.createElement("p");
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
let renderer;
let system = new System("https://arolauntech.github.io/kspmga/data/systems/stock.json", (sys, json) => {
	renderer = new Renderer({
		fov: 75,
		scale: json.SCALE,
		mindist: 0.01,
		maxdist: 100000,
		aspect: aspectratio,
		minscrolldist: json.mindist,
		maxscrolldist: json.maxdist,
		htmlCanvas: canvas
	});

	renderer.fillSceneWithSystem(sys);
	renderer.updateSceneWithSystem(sys, 0);

	setBodySelectOptions(initbody, sys);
	setBodySelectOptions(finalbody, sys);

	setTimeSystem(json.timesystem);
});

let mgafinder;