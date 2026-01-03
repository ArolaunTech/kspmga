import { Orbit } from './orbit.js';
import { Vector3 } from './vector.js';

export class System {
	bodies = [];
	bodymap = new Map();
	abbreviations = new Map();

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

				this.generateAbbreviations();

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

	generateAbbreviations() {
		let massindices = [];
		for (let i = 0; i < this.bodies.length; i++) {
			massindices.push(i);
		}

		let bodies = this.bodies;
		massindices.sort(function (a, b) {
			return bodies[b].gravparameter - bodies[a].gravparameter;
		});

		let maxlength = 0;
		for (let i = 0; i < this.bodies.length; i++) {
			if (this.bodies[i].name.length > maxlength) {
				maxlength = this.bodies[i].name.length;
			}
		}

		for (let i = 0; i < this.bodies.length; i++) {
			this.abbreviations.set(this.bodies[i].name, this.bodies[i].name);
		}
		for (let length = 1; length <= maxlength; length++) {
			for (let i = 0; i < massindices.length; i++) {
				if (length >= this.bodies[massindices[i]].name.length) continue;

				let abbreviation = this.bodies[massindices[i]].name.slice(0, length);

				if (this.abbreviations.has(abbreviation)) continue;

				this.abbreviations.set(abbreviation, this.bodies[massindices[i]].name);
			}
		}
	}

	getSequenceRegex() {
		return "^((" + 
			Array.from(this.abbreviations.keys()).join("|") +
			")-)+(" +
			Array.from(this.abbreviations.keys()).join("|") +
			")$";
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
};