// Elements
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
const timeslider = document.getElementById("timeslider");
const timedisplay = document.getElementById("timelabel");

// Functions
export function handleNumericMinMax() {
	if (this.value === "") return;

	let numericvalue = Number(this.value);

	if ((this.min !== "") && (numericvalue < Number(this.min))) {
		this.value = this.min;
	}

	if ((this.max !== "") && (numericvalue > Number(this.max))) {
		this.value = this.max;
	}
}

export function handleNumericMinMaxInt() {
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

export function setBodySelectOptions(element, sys) {
	element.textContent = "";

	for (let i = 0; i < sys.bodies.length; i++) {
		if (sys.bodies[i].root) continue; // Cannot currently do flyby sequences with arbitrary solar orbits

		const option = document.createElement("option");

		option.innerText = sys.bodies[i].name;

		element.appendChild(option);
	}
}

// Actions
[initalt, finalalt, minvinf, maxvinf, maxdvdsm, maxduration].forEach((element) => {
	element.onchange = handleNumericMinMax;
});

[earliestyear, earliesthour, earliestday, latestyear, latesthour, latestday].forEach((element) => {
	element.onchange = handleNumericMinMaxInt;
});