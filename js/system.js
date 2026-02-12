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
				//console.log(structuredClone(json));

				this.bodies = json.bodies;

				let numbodies = this.bodies.length;
				for (let i = 0; i < numbodies; i++) {
					this.bodymap.set(this.bodies[i].name, i);
				}

				for (let i = 0; i < numbodies; i++) {
					if (this.bodies[i].root) continue;

					//console.log(i, this.bodies[i]);

					let parent = this.bodies[i].orbit.parent;
					let mu = this.bodies[this.bodymap.get(parent)].gravparameter;

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
};