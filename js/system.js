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

				this.ready = true;
			});
	}

	getPosition(name, time) {
		let body = this.bodies[this.bodymap.get(name)];

		if (body.root) {
			
		}
	}
};