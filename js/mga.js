import { System } from './system.js';
import { Vector3 } from './vector.js';
import { realRootsQuartic } from './polynomial.js';
import { propagate } from './orbit.js';
import { TrajectoryTree } from './tree.js';
import { randint } from './random.js';

//Compute optimal MGA trajectories
function infmult(M, ga) {
	if (M === 0) return 0;
	return M * ga;
}

export class MGAFinder {
	constructor(system) {
		this.system = system;
	}

	// Find initial velocities in the compatibility method.
	findInitialVelocities(r1, v1, r2, rv, mu) {
		let dist1 = r1.norm;
		let dist2 = r2.norm;
		let rdot = Vector3.dot(r1, r2);
		let rdiff = Vector3.sub(r2, r1);
		let c = rdiff.norm;
		let uc = Vector3.normalize(rdiff);
		let ui = Vector3.normalize(r1);

		let Kh = mu * c / (dist1 * dist2 + rdot);

		if ((Kh === 0) || ((dist1 * dist2 + rdot) === 0)) {
			// n-pi transfer - there are an infinite number of possible velocities.
			// This is handled elsewhere.
			return [];
		}

		let rootKh = Math.sqrt(Kh);

		let P = -2 * Vector3.dot(ui, v1) / rootKh;
		let Q = Vector3.dot(v1, v1) / Kh - rv * rv / Kh + 2 * Vector3.dot(ui, uc);
		let R = -2 * Vector3.dot(uc, v1) / rootKh;
		let S = 1;

		let vps = realRootsQuartic(P, Q, R, S);

		let vels = [];

		for (let i = 0; i < vps.length; i++) {
			if (vps[i] === 0) {
				continue;
			}

			vels.push(Vector3.add(
				Vector3.mult(ui, rootKh * vps[i]),
				Vector3.mult(uc, rootKh / vps[i])
			));
		}

		return vels;
	}

	Fa(r1, vout, r2, t1, t2, mu) {
		let rdiff = Vector3.sub(r2, r1);
		let uc = Vector3.normalize(rdiff);
		let ui = Vector3.normalize(r1);
		let u = Vector3.add(ui, uc);

		let dist1 = r1.norm;
		let dist2 = r2.norm;
		let vel1 = vout.norm;

		let energy = vel1 * vel1 / 2 - mu / dist1;
		let sma = -mu / 2 / energy;

		let h = Vector3.cross(r1, vout).norm;
		let p = h * h / mu;

		let rvdot = Vector3.dot(r1, vout);

		let cosdE = 1 - (dist1 * dist2 - Vector3.dot(r1, r2)) / (sma * p);
		let sindE = Math.sqrt(Math.abs(1 - cosdE * cosdE));

		let dE;

		let cosdTheta = Vector3.dot(r1, r2) / dist1 / dist2;

		if (energy < 0) {
			// Elliptical case

			//   cos(dEa)
			// = 1 - (1 - cos(dEa))
			// = 1 - (1 - cos(dEa))(1 - e^2) / (1 - e^2)
			// = 1 - ((1 - e^2) - (1 - e^2)cos(dEa)) / (1 - e^2)
			// = 1 - ((1 - e^2) + (e^2 - 1)cos(dEa)) / (1 - e^2)
			// = 1 - (1 - e^2 + e^2cos(dEa) - cos(dEa)) / (1 - e^2)
			// = 1 - (1 - e^2 + e^2cos(E0)cos(E1) + e^2sin(E0)sin(E1) - cos(E0)cos(E1) - sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (1 - e^2 - ecos(E0) - ecos(E1) + e^2cos(E0)cos(E1) + ecos(E0) + ecos(E1) + e^2sin(E0)sin(E1) - cos(E0)cos(E1) - sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - ((1 - ecos(E0)) - ecos(E1)(1 - ecos(E0)) - e^2 + ecos(E0) + ecos(E1) + e^2sin(E0)sin(E1) - cos(E0)cos(E1) - sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - ((1 - ecos(E0))(1 - ecos(E1)) - e^2 + ecos(E0) + ecos(E1) + e^2sin(E0)sin(E1) - cos(E0)cos(E1) - sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - e^2 + ecos(E0) + ecos(E1) + e^2sin(E0)sin(E1) - cos(E0)cos(E1) - sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - e^2 + ecos(E0) + ecos(E1) - cos(E0)cos(E1) - (1 - e^2)sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - (e^2 - ecos(E0) - ecos(E1) + cos(E0)cos(E1)) - (1 - e^2)sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - (e^2 - ecos(E1) + cos(E0)cos(E1) - ecos(E0)) - (1 - e^2)sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - (e(e - cos(E1)) + cos(E0)(cos(E1) - e)) - (1 - e^2)sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - (cos(E0) - e)(cos(E1) - e) - (1 - e^2)sin(E0)sin(E1)) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - (cos(E0) - e)(cos(E1) - e) - b^2 * sin(E0)sin(E1) / a^2) / (1 - e^2)
			// = 1 - (r0r1 / a^2 - a * (cos(E0) - e) * a * (cos(E1) - e) / a^2 - b^2 * sin(E0)sin(E1) / a^2) / (1 - e^2)
			// = 1 - (r0r1 - a * (cos(E0) - e) * a * (cos(E1) - e) - b^2 * sin(E0)sin(E1)) / (a^2 * (1 - e^2))
			// = 1 - (r0r1 - a * (cos(E0) - e) * a * (cos(E1) - e) - b^2 * sin(E0)sin(E1)) / (ap)
			// = 1 - (r0r1 - x0 * x1 - y0 * y1) / (ap)
			// 
			// +================ Finding dEa ================+
			// | cos(dEa) = 1 - (r0*r1 - dot(r0, r1)) / (ap) |
			// +=============================================+

			if ((dist1 * dist2 - Vector3.dot(r1, r2)) > (dist1 + dist2) * p) {
				sindE *= -1;
			}

			if (Vector3.dot(vout, u) < 0) {
				sindE *= -1;
			}

			dE = Math.atan2(sindE, cosdE);
			if (dE < 0) {
				dE += 2 * Math.PI;
			}
		} else {
			// Hyperbolic case

			// cosh(dF)
			// = 1 - (1 - cosh(dF))
			// = 1 - (1 - cosh(dF))(e^2 - 1) / (e^2 - 1)
			// = 1 - (e^2 - 1 - (e^2 - 1)cosh(dF)) / (e^2 - 1)
			// = 1 - (e^2 - 1 - e^2 * cosh(dF) + cosh(dF)) / (e^2 - 1)
			// = 1 - (e^2 - 1 - e^2 * cosh(F0)cosh(F1) + e^2sinh(F0)sinh(F1) + cosh(F0)cosh(F1) - sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - (e^2 - 1 - e^2 * cosh(F0)cosh(F1) + cosh(F0)cosh(F1) + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - (e^2 - 1 - ecosh(F0) - ecosh(F1) - e^2 * cosh(F0)cosh(F1) + ecosh(F0) + ecosh(F1) + cosh(F0)cosh(F1) + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - (e^2 - ecosh(F0) - ecosh(F1) + cosh(F0)cosh(F1) - 1 + ecosh(F0) + ecosh(F1) - e^2 * cosh(F0)cosh(F1) + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - ((e - cosh(F0))(e - cosh(F1)) - 1 + ecosh(F0) + ecosh(F1) - e^2 * cosh(F0)cosh(F1) + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - ((e - cosh(F0))(e - cosh(F1)) - (1 - ecosh(F0))(1 - ecosh(F1)) + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - ((e - cosh(F0))(e - cosh(F1)) - r0 * r1 / a^2 + (e^2 - 1)sinh(F0)sinh(F1)) / (e^2 - 1)
			// = 1 - ((e - cosh(F0))(e - cosh(F1)) - r0 * r1 / a^2 + b^2*sinh(F0)sinh(F1)/a^2) / (e^2 - 1)
			// = 1 - (a * (e - cosh(F0)) * a * (e - cosh(F1)) - r0 * r1 + b^2*sinh(F0)sinh(F1)) / (a^2 * (e^2 - 1))
			// = 1 - (a * (cosh(F0) - e) * a * (cosh(F1) - e) - r0 * r1 + b^2*sinh(F0)sinh(F1)) / (a^2 * (e^2 - 1))
			// = 1 - (x0 * x1 - r0 * r1 + y0 * y1) / (a^2 * (e^2 - 1))
			// = 1 - (dot(r0, r1) - r0 * r1) / (a^2 * (e^2 - 1))
			// = 1 - (r0 * r1 - dot(r0, r1)) / (a^2 * (1 - e^2))
			// 
			// +================ Finding dF =================+
			// | cosh(dF) = 1 - (r0*r1 - dot(r0, r1)) / (ap) |
			// +=============================================+

			// Because this formula is the exact same as for the elliptic case, cosh(dF)
			// will be in the variable cosdE.

			// TO-DO: check if code is correct.
			let dTheta = Math.acos(cosdTheta);
			if (Vector3.dot(vout, u) < 0) {
				dTheta = 2 * Math.PI - dTheta;
			}

			let eccentricity = Math.sqrt(1 - p / sma);
			let f = Math.acos((p / dist1 - 1) / eccentricity);

			if (rvdot < 0) {
				f *= -1;
			}

			let maxdTheta = Math.acos(-1 / eccentricity) - f;
			if (dTheta > maxdTheta) {
				return -Infinity;
			}

			dE = Math.asinh(sindE);
		}

		// This formula is slightly different than the one in Gavira-Aladro and Bombardelli (2024), but it gives the same answer
		let dt = sma * Math.sqrt(Math.abs(sma) / mu) * (rvdot * (1 - cosdE) / Math.sqrt(mu * Math.abs(sma)) - (1 - dist1 / sma) * sindE + dE);

		return t2 - t1 - dt;
	}

	Ga(r1, vout, mu) {
		let dist1 = r1.norm;
		let vel = vout.norm;

		let energy = vel * vel / 2 - mu / dist1;

		if (energy >= 0) {
			return Infinity;
		}

		let sma = -mu / 2 / energy;

		return 2 * Math.PI * sma * Math.sqrt(sma / mu);
	}

	findTransfersInt(r1, v1, orbit, t1, rv, mu, maxrevs) {
		const discretizationSteps = 100;
		const epsilon = 1e-3;

		let period = orbit.period;

		let fas = [[], [], [], []];
		let gas = [[], [], [], []];

		let found = 0;

		let minfas = [Infinity, Infinity, Infinity, Infinity];
		let mingas = [Infinity, Infinity, Infinity, Infinity];
		let maxfas = [-Infinity, -Infinity, -Infinity, -Infinity];
		let maxgas = [-Infinity, -Infinity, -Infinity, -Infinity];

		for (let i = 0; i < discretizationSteps; i++) {
			let t2 = t1 + period * i / discretizationSteps;

			let state2 = orbit.getState(t2);
			let vels = this.findInitialVelocities(r1, v1, state2.r, rv, mu);

			if (vels.length > found) {
				found = vels.length;
			}

			for (let j = 0; j < vels.length; j++) {
				let fa = this.Fa(r1, vels[j], state2.r, t1, t2, mu);
				let ga = this.Ga(r1, vels[j], mu);

				if ((fa > maxfas[j]) && Number.isFinite(fa)) maxfas[j] = fa;
				if (ga > maxgas[j]) maxgas[j] = ga;
				if ((fa < minfas[j]) && Number.isFinite(fa)) minfas[j] = fa;
				if (ga < mingas[j]) mingas[j] = ga;

				fas[j].push(fa);
				gas[j].push(ga);
			}

			for (let j = vels.length; j < 4; j++) {
				fas[j].push(null);
				gas[j].push(null);
			}
		}

		let out = [];

		for (let i = 0; i < found; i++) {
			for (let N = 0; N < maxrevs; N++) {
				for (let M = 0; M < maxrevs; M++) {
					let mincompat = minfas[i] + N * period - infmult(M, maxgas[i]);
					let maxcompat = maxfas[i] + N * period - infmult(M, mingas[i]);

					if (mincompat > 0) continue;
					if (maxcompat < 0) break;

					for (let j = 0; j < discretizationSteps - 1; j++) {
						if (!Number.isFinite(fas[i][j])) continue;
						if (!Number.isFinite(fas[i][j + 1])) continue;

						if (gas[i][j] === null) continue;
						if (gas[i][j + 1] === null) continue;

						if ((M > 0) && !(Number.isFinite(gas[i][j]) && Number.isFinite(gas[i][j + 1]))) continue;

						let compat = fas[i][j] + N * period - infmult(M, gas[i][j]);
						let compatnext = fas[i][j + 1] + N * period - infmult(M, gas[i][j + 1]);

						if (compat * compatnext > 0) continue;

						let lowx = t1 + period * (j / discretizationSteps + N);
						let lowy = compat;

						let highx = t1 + period * ((j + 1) / discretizationSteps + N);
						let highy = compatnext;

						let mid;
						let midy = 2;

						let vels;

						// Regula falsi
						for (let k = 0; (k < 10) && Math.abs(midy) >= 1; k++) {
							let slope = (highy - lowy) / (highx - lowx);
							mid = lowx - lowy / slope;

							let state2 = orbit.getState(mid);
							vels = this.findInitialVelocities(r1, v1, state2.r, rv, mu);

							if (i >= vels.length) break;

							let fa = this.Fa(r1, vels[i], state2.r, t1, mid, mu);
							let ga = this.Ga(r1, vels[i], mu);

							midy = fa - infmult(M, ga);

							if (midy * lowy > 0) {
								lowx = mid;
								lowy = midy;
							} else {
								highx = mid;
								highy = midy;
							}
						}

						if (i >= vels.length) continue;

						//console.log(mid, midy);

						let state2 = orbit.getState(mid);

						let outgoing = Vector3.sub(vels[i], v1);
						let incoming = Vector3.sub(propagate(r1, vels[i], mid - t1, mu).v, state2.v);

						out.push({
							t1: t1,
							t2: mid,
							outgoing: outgoing,
							outgoingmag: outgoing.norm,
							incoming: incoming,
							incomingmag: incoming.norm
						});
					}
				}
			}
		}

		return out;
	}

	// Find no-DSM transfers betweem two planets.
	findTransfersNoDSM(planet1, planet2, t1, rv, maxrevs) {
		let state1 = this.system.getDState(planet1, t1);

		let parent = this.system.bodies[this.system.bodymap.get(planet1)].parent;

		return this.findTransfersInt(state1.r, state1.v, this.system.bodies[this.system.bodymap.get(planet2)].orbit, t1, rv, this.system.bodies[this.system.bodymap.get(parent)].gravparameter, maxrevs);
	}

	calcEjectiondV(body, vinf, circalt) {
		let stdgp = this.system.bodies[this.system.bodymap.get(body)].gravparameter;
		let soi = this.system.bodies[this.system.bodymap.get(body)].soi;

		let energy = 0.5 * vinf * vinf - stdgp / soi;
		let sma = -stdgp / 2 / energy;

		let ejectionvel = Math.sqrt(stdgp * (2 / circalt - 1 / sma));
		let circvel = Math.sqrt(stdgp / circalt);

		return ejectionvel - circvel;
	}

	calcCosDeflection(body, vinf) {
		let stdgp = this.system.bodies[this.system.bodymap.get(body)].gravparameter;
		let minalt = this.system.bodies[this.system.bodymap.get(body)].radius + this.system.bodies[this.system.bodymap.get(body)].atmodepth + 10000;

		return 1 - 2 * Math.pow(1 / (1 + minalt * vinf * vinf / stdgp), 2);
	}

	planMGATrajectory(sequence, initalt, finalalt, minvinf, maxvinf, maxdvdsm, maxduration, maxrevs, earliesttime, latesttime, includecapture) {
		let fullsequence = [];
		for (let i = 0; i < sequence.length; i++) {
			fullsequence.push(this.system.abbreviations.get(sequence[i]));
		}

		let trueinitalt = this.system.bodies[this.system.bodymap.get(fullsequence[0])].radius + initalt;
		let truefinalalt = this.system.bodies[this.system.bodymap.get(fullsequence[fullsequence.length - 1])].radius + finalalt;

		let truelatesttime = latesttime;
		if (latesttime === Infinity) {
			truelatesttime = earliesttime + 1e9; // About 10 years
		}

		console.log(fullsequence, trueinitalt, truefinalalt, minvinf, maxvinf, maxdvdsm, maxduration, maxrevs, earliesttime, truelatesttime, includecapture);

		let tree = new TrajectoryTree(0, 0, {});

		console.log(tree);

		for (let steps = 0; steps < 100; steps++) {
			let randomnode = tree;

			while (randomnode.children.length > 0) {
				if (Math.random < 0.2) break;
				if (randomnode.planet === fullsequence.length - 2) break;

				randomnode = randomnode.children[randint(0, randomnode.children.length)];
			}

			if (randomnode.planet === 0) {
				let vinf = minvinf + Math.random() * (maxvinf - minvinf);
				let start = earliesttime + Math.random() * (truelatesttime - earliesttime);

				let transfers = this.findTransfersNoDSM(fullsequence[0], fullsequence[1], start, vinf, maxrevs);

				for (let i = 0; i < transfers.length; i++) {
					randomnode.children.push(new TrajectoryTree(1, this.calcEjectiondV(fullsequence[0], transfers[i].outgoingmag, trueinitalt), transfers[i]));
				}
			} else {
				let vinf = randomnode.prevtransfer.incomingmag;
				let start = randomnode.prevtransfer.t2;

				let transfers = this.findTransfersNoDSM(fullsequence[0], fullsequence[1], start, vinf, maxrevs);

				for (let i = 0; i < transfers.length; i++) {
					let cosdeflection = Vector3.dot(transfers[i].outgoing, randomnode.prevtransfer.incoming) / vinf / vinf;
					let mincosdeflection = this.calcCosDeflection(fullsequence[randomnode.planet], vinf);

					if (cosdeflection < mincosdeflection) continue;

					randomnode.children.push(new TrajectoryTree(randomnode.planet + 1, randomnode.dvused, transfers[i]));
				}
			}
		}

		console.log(tree);

		let paths = tree.getAllLeafPaths();

		let trajectories = [];

		let mindv = Infinity;

		for (let i = 0; i < paths.length; i++) {
			if (paths[i].length === fullsequence.length) {
				paths[i].push(
					paths[i][fullsequence.length - 1].dvused + 
					this.calcEjectiondV(
						fullsequence[fullsequence.length - 1],
						paths[i][fullsequence.length - 1].prevtransfer.incomingmag,
						truefinalalt
					)
				);

				if (paths[i][fullsequence.length] < mindv) {
					console.log(paths[i]);

					mindv = paths[i][fullsequence.length];
				}

				trajectories.push(paths[i]);
			}
		}

		console.log(trajectories);
	}
};