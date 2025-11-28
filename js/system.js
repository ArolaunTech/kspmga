class System {
	bodies = [];
	bodymap = new Map();

	constructor(JSONPath) {
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

				this.ready = true;
			});
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
	}
};