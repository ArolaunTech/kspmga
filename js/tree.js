export class TrajectoryTree {
	constructor(planet, dvused, prevtransfer) {
		this.children = [];
		this.planet = planet;
		this.dvused = dvused;
		this.prevtransfer = prevtransfer;
	}

	getAllLeafPaths() {
		if (this.children.length === 0) {
			return [[this]];
		}

		let out = [];

		for (let i = 0; i < this.children.length; i++) {
			let paths = this.children[i].getAllLeafPaths();

			for (let j = 0; j < paths.length; j++) {
				out.push([this].concat(paths[j]));
			}
		}

		return out;
	}
};