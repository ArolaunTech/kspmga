import { System } from './system.js';
import { Vector3 } from './vector.js';
import { realRootsQuartic } from './polynomial.js';
import { Orbit, propagate } from './orbit.js';
import { TrajectoryTree } from './tree.js';
import { randint, normal } from './random.js';
import { lambert } from './lambert.js';

//Compute optimal MGA trajectories
function infmult(M, ga) {
	if (M === 0) return 0;
	return M * ga;
}

class MGAFinder {
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
			for (let M = 0; M <= maxrevs; M++) {
				for (let N = 0; N <= 100; N++) {
					let mincompat = minfas[i] + N * period - infmult(M, maxgas[i]);
					let maxcompat = maxfas[i] + N * period - infmult(M, mingas[i]);

					if (mincompat > 0) break;
					if (maxcompat < 0) continue;

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

							if (slope === 0) break;

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

						if (!Number.isFinite(mid)) continue;
						if (Vector3.sub(propagate(r1, vels[i], mid - t1, mu).r, state2.r).norm > 1e4) continue;
						if (!Number.isFinite(outgoing.norm)) continue; // Resonant orbits can cause this - I saw this cause an error with a transfer that encountered Moho after 6 revolutions of Moho
						if (!Number.isFinite(incoming.norm)) continue;

						//if (Math.abs((mid - t1) / 9201600 - 0.5) < 0.1)
						//	console.log(Vector3.sub(propagate(r1, vels[i], mid - t1, mu).r, state2.r).norm, mid - t1, (mid - t1) / 9201600, outgoing.norm, incoming.norm);

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

	violation(body, vin, vout) {
		let innorm = vin.norm;
		let outnorm = vout.norm;

		let cosmaxdeflection = this.calcCosDeflection(body, innorm);
		let cosdeflection = Vector3.dot(vin, vout) / innorm / outnorm;

		if (cosdeflection > cosmaxdeflection) {
			return Math.abs(outnorm - innorm);
		}

		let sindir = Vector3.normalize(Vector3.cross(Vector3.cross(vin, vout), vin));
		let closest = Vector3.add(Vector3.mult(vin, cosmaxdeflection), Vector3.mult(sindir, innorm * Math.sqrt(1 - cosmaxdeflection * cosmaxdeflection)));

		return Vector3.sub(closest, vout).norm;
	}

	mutate(dir, jiggle) {
		let out = Vector3.add(dir, new Vector3(jiggle * normal(), jiggle * normal(), jiggle * normal()));

		return Vector3.normalize(out);
	}

	fitness(trajectory, sequence, initalt, finalalt, includecapture) {
		let parent = this.system.bodies[this.system.bodymap.get(sequence[0])].parent;
		let mu = this.system.bodies[this.system.bodymap.get(parent)].gravparameter;

		let outgoing = new Vector3(trajectory[1], trajectory[2], trajectory[3]);

		let dv = this.calcEjectiondV(sequence[0], outgoing.norm, initalt);

		// Iterate over each leg
		for (let i = 0; i < sequence.length - 1; i++) {
			let state1 = this.system.getDState(sequence[i], trajectory[5 * i]);
			let state2 = this.system.getDState(sequence[i + 1], trajectory[5 * i + 5]);
			let totaltime = trajectory[5 * i + 5] - trajectory[5 * i];

			let vout = Vector3.add(
				state1.v, 
				new Vector3(
					trajectory[5 * i + 1],
					trajectory[5 * i + 2],
					trajectory[5 * i + 3]
				));

			let statedsm = propagate(state1.r, vout, trajectory[5 * i + 4] * totaltime, mu);

			let sols = lambert(statedsm.r, state2.r, (1 - trajectory[5 * i + 4]) * totaltime, mu, 5, false);

			let mindsmdv = Infinity;
			let bestvin;

			for (let j = 0; j < sols.length; j++) {
				let vafter = sols[j][0];
				let dsmdv = Vector3.sub(vafter, statedsm.v).norm;

				if (dsmdv < mindsmdv) {
					mindsmdv = dsmdv;
					bestvin = sols[j][1];
				}
			}

			dv += mindsmdv;

			let bestincoming = Vector3.sub(bestvin, state2.v);

			if (i === sequence.length - 2) {
				// Capture
				if (includecapture) {
					dv += this.calcEjectiondV(sequence[i + 1], bestincoming.norm, finalalt);
				}
			} else {
				// Flyby
				let postflybyoutgoing = new Vector3(
					trajectory[5 * i + 6],
					trajectory[5 * i + 7],
					trajectory[5 * i + 8]
				);

				dv += this.violation(sequence[i + 1], bestincoming, postflybyoutgoing);
			}
		}

		return dv;
	}

	planMGATrajectory(sequence, initalt, finalalt, minvinf, maxvinf, maxdvdsm, maxduration, maxrevs, earliesttime, latesttime, includecapture, flybydvs) {
		let trueinitalt = this.system.bodies[this.system.bodymap.get(sequence[0])].radius + initalt;
		let truefinalalt = this.system.bodies[this.system.bodymap.get(sequence[sequence.length - 1])].radius + finalalt;

		let truelatesttime = latesttime;
		if (latesttime === Infinity) {
			truelatesttime = earliesttime + 1e9; // About 10 years
		}

		let tree = new TrajectoryTree(0, 0, {});
		let jiggle = maxdvdsm / 5000;
		let timejiggle = 5000 * maxdvdsm;

		let parent = this.system.bodies[this.system.bodymap.get(sequence[0])].parent;
		let mu = this.system.bodies[this.system.bodymap.get(parent)].gravparameter;

		// Global search
		let numnodes = 0;
		let totalnodes = 0;
		for (let steps = 0; steps < 50000; steps++) {
			if (steps % 1000 === 0) {
				postMessage({
					status: "report",
					steps: steps,
					foundtrajectories: numnodes,
					nodes: totalnodes
				});
			}
			
			if (numnodes > 2500) break;
			if ((numnodes > 10) && (totalnodes > 10000)) break;

			let randomnode = tree;

			while (randomnode.children.length > 0) {
				if (Math.random() < 1 / 5) break;
				if (randomnode.planet === sequence.length - 2) break;

				randomnode = randomnode.children[randint(0, randomnode.children.length)];
			}

			if (randomnode.children.length > steps / 1000 && ((randomnode.planet > 0) || (Math.random() < 0.9))) {
				steps--;
				continue;
			}

			//console.log(randomnode.planet);

			let resonant = sequence[randomnode.planet] === sequence[randomnode.planet + 1];

			if (randomnode.planet === 0) {
				if (Math.random() < 0.99) steps--;

				let vinf = minvinf + Math.random() * (maxvinf - minvinf);
				let start = earliesttime + Math.random() * (truelatesttime - earliesttime);

				//console.log(vinf, start);

				let transfers = this.findTransfersNoDSM(sequence[0], sequence[1], start, vinf, maxrevs[0]);

				for (let i = 0; i < transfers.length; i++) {
					//if (Math.random() < 0.9) continue;

					randomnode.children.push(new TrajectoryTree(1, this.calcEjectiondV(sequence[0], vinf, trueinitalt), transfers[i]));
				
					if (sequence.length === 2) numnodes++;
					totalnodes++;
				}

				if (resonant) {
					let state1 = this.system.getDState(sequence[0], start);
					let vel = state1.v.norm;
					let dist1 = state1.r.norm;
	
					let maxvel = vel + vinf;
					let minvel = Math.max(0, vel - vinf);
	
					let maxsma = 1 / (2 / dist1 - maxvel * maxvel / mu);
					let minsma = 1 / (2 / dist1 - maxvel * maxvel / mu);
	
					let maxperiod;
	
					if (maxsma <= 0) {
						maxperiod = Infinity;
					} else {
						maxperiod = 2 * Math.PI * maxsma * Math.sqrt(maxsma / mu);
					}
					let minperiod = 2 * Math.PI * minsma * Math.sqrt(minsma / mu);
	
					let currperiod = this.system.bodies[this.system.bodymap.get(sequence[0])].orbit.period;
	
					for (let b = 1; b <= maxrevs[0]; b++) {
						for (let a = 1; a <= 100; a++) {
							let period = currperiod * a / b;
	
							if (period > maxperiod) break;
							if (period < minperiod) continue;
	
							let sma = Math.cbrt(mu * Math.pow(period / 2 / Math.PI, 2));
							let targetvel = Math.sqrt(mu * (2 / dist1 - 1 / sma));
	
							let vforward = (targetvel * targetvel - vinf * vinf - vel * vel) / (2 * vel);
							let vside = Math.sqrt(vinf * vinf - vforward * vforward);

							if (!Number.isFinite(vside)) continue;
	
							let forwarddir = Vector3.normalize(state1.v);
							let sidedir1 = Vector3.normalize(Vector3.cross(new Vector3(0, 0, 1), forwarddir));
							let sidedir2 = Vector3.normalize(Vector3.cross(sidedir1, forwarddir));
	
							let angle = Math.random() * 2 * Math.PI;
	
							let outgoing = Vector3.add(
								Vector3.mult(forwarddir, vforward), 
								Vector3.add(
									Vector3.mult(sidedir1, vside * Math.cos(angle)), 
									Vector3.mult(sidedir2, vside * Math.sin(angle))
								)
							);
	
							randomnode.children.push(
								new TrajectoryTree(
									1, 
									this.calcEjectiondV(sequence[0], vinf, trueinitalt), 
									{
										t1: start,
										t2: start + currperiod * a, 
										outgoing: outgoing, 
										outgoingmag: vinf, 
										incoming: outgoing, 
										incomingmag: vinf
									}
								)
							);

							if (sequence.length === 2) numnodes++;
							totalnodes++;
						}
					}
				}
			} else {
				let vinf = randomnode.prevtransfer.incomingmag;
				let maxflybydv = flybydvs[randomnode.planet - 1] + 1;
				vinf += normal() * maxflybydv * 0.5;

				let previncoming = randomnode.prevtransfer.incoming;
				let start = randomnode.prevtransfer.t2;

				let transfers = this.findTransfersNoDSM(sequence[randomnode.planet], sequence[randomnode.planet + 1], start, vinf, maxrevs[randomnode.planet]);

				for (let i = 0; i < transfers.length; i++) {
					//if (Math.random() < 0.9) continue;
					let nonjiggledviolation = this.violation(sequence[randomnode.planet], previncoming, transfers[i].outgoing);
					if (nonjiggledviolation < maxflybydv) {
						randomnode.children.push(new TrajectoryTree(randomnode.planet + 1, randomnode.dvused + nonjiggledviolation, transfers[i]));
					
						if (randomnode.planet === sequence.length - 2) numnodes++;
						totalnodes++;
					}
			
					if (maxdvdsm > 0) {
						let nodsmoutdir = Vector3.normalize(transfers[i].outgoing);
						let jiggled = Vector3.mult(this.mutate(nodsmoutdir, jiggle), vinf);
						let dvdsm = Math.random() * maxdvdsm;

						let propagatetime = Math.random() * (transfers[i].t2 - transfers[i].t1);

						let state1 = this.system.getDState(sequence[randomnode.planet], start);
						let stateDSM = propagate(state1.r, Vector3.add(state1.v, jiggled), propagatetime, mu);

						let jiggledviolation = this.violation(sequence[randomnode.planet], previncoming, jiggled);
						if (jiggledviolation > maxflybydv) continue;

						let t2 = transfers[i].t2 + normal() * timejiggle;

						let state2 = this.system.getDState(sequence[randomnode.planet + 1], t2);

						//console.log(transfers[i]);

						let sols = lambert(stateDSM.r, state2.r, t2 - start - propagatetime, mu, 5, false);

						let mindvdsm = Infinity;
						let bestvout;
						let bestvin;

						for (let j = 0; j < sols.length; j++) {
							let dvdsm = Vector3.sub(sols[j][0], stateDSM.v).norm;

							if (dvdsm < mindvdsm) {
								mindvdsm = dvdsm;
								bestvout = sols[j][0];
								bestvin = sols[j][1];
							}
						}

						if (mindvdsm > maxdvdsm) continue;

						let dsmtransfer = {
							t1: start,
							t2: t2,
							tdsm: start + propagatetime,
							rdsm: stateDSM.r,
							outgoing: jiggled,
							outgoingmag: vinf,
							dsm: Vector3.sub(bestvout, stateDSM.v),
							dsmmag: mindvdsm,
							incoming: Vector3.sub(bestvin, state2.v),
							incomingmag: Vector3.sub(bestvin, state2.v).norm
						};

						randomnode.children.push(new TrajectoryTree(randomnode.planet + 1, randomnode.dvused + mindvdsm + jiggledviolation, dsmtransfer));
								
						if (randomnode.planet === sequence.length - 2) numnodes++;
						totalnodes++;
					}
				}

				if (resonant) {
					let state1 = this.system.getDState(sequence[randomnode.planet], start);
					let vel = state1.v.norm;
					let dist1 = state1.r.norm;

					let minvel = Math.max(0, vel - vinf);
					let minsma = 1 / (2 / dist1 - minvel * minvel / mu);
					let minperiod = 2 * Math.PI * minsma * Math.sqrt(minsma / mu);
	
					let currperiod = this.system.bodies[this.system.bodymap.get(sequence[randomnode.planet])].orbit.period;
	
					for (let b = 1; b <= maxrevs[randomnode.planet]; b++) {
						for (let a = 1; a <= 100; a++) {
							let period = currperiod * a / b;
	
							if (period < minperiod) continue;
	
							let sma = Math.cbrt(mu * Math.pow(period / 2 / Math.PI, 2));

							if (sma < dist1 / 2) continue; // Periapsis < 0

							let targetvel = Math.sqrt(mu * (2 / dist1 - 1 / sma));
	
							let vforward = (targetvel * targetvel - vinf * vinf - vel * vel) / (2 * vel);

							if (vforward > vinf) break;

							let vside = Math.sqrt(vinf * vinf - vforward * vforward);

							if (!Number.isFinite(vside)) continue;
	
							let forwarddir = Vector3.normalize(state1.v);
							let sidedir1 = Vector3.normalize(Vector3.cross(new Vector3(0, 0, 1), forwarddir));
							let sidedir2 = Vector3.normalize(Vector3.cross(sidedir1, forwarddir));
	
							let angle = Math.random() * 2 * Math.PI;
	
							let nodsmoutdir = Vector3.add(
								Vector3.mult(forwarddir, vforward / vinf), 
								Vector3.add(
									Vector3.mult(sidedir1, vside / vinf * Math.cos(angle)), 
									Vector3.mult(sidedir2, vside / vinf * Math.sin(angle))
								)
							);

							let jiggled = Vector3.mult(this.mutate(nodsmoutdir, jiggle), vinf);

							let jiggledviolation = this.violation(sequence[randomnode.planet], previncoming, jiggled);
							if (jiggledviolation > maxflybydv) continue;
	
							let propagatetime = Math.random() * currperiod * a;

							let stateDSM = propagate(state1.r, Vector3.add(state1.v, jiggled), propagatetime, mu);
	
							let t2 = start + currperiod * a;
							let tdsm = start + propagatetime;

							let t2jiggled = start + currperiod * a + normal() * timejiggle;

							let state2 = this.system.getDState(sequence[randomnode.planet + 1], t2jiggled);

							let sols = lambert(stateDSM.r, state2.r, t2jiggled - tdsm, mu, 5, false);

							let mindvdsm = Infinity;
							let bestvout;
							let bestvin;

							for (let j = 0; j < sols.length; j++) {
								let dvdsm = Vector3.sub(sols[j][0], stateDSM.v).norm;

								if (dvdsm < mindvdsm) {
									mindvdsm = dvdsm;
									bestvout = sols[j][0];
									bestvin = sols[j][1];
								}
							}

							if (mindvdsm > maxdvdsm) continue;

							let dsmtransfer = {
								t1: start,
								t2: t2jiggled,
								tdsm: tdsm,
								rdsm: stateDSM.r,
								outgoing: jiggled,
								outgoingmag: vinf,
								dsm: Vector3.sub(bestvout, stateDSM.v),
								dsmmag: mindvdsm,
								incoming: Vector3.sub(bestvin, state2.v),
								incomingmag: Vector3.sub(bestvin, state2.v).norm
							};
	
							randomnode.children.push(new TrajectoryTree(randomnode.planet + 1, randomnode.dvused + mindvdsm + jiggledviolation, dsmtransfer));
							
							if (randomnode.planet === sequence.length - 2) numnodes++;
							totalnodes++;
						}
					}
				}
			}
		}

		// Get all found trajectories
		let paths = tree.getAllLeafPaths();

		let trajectories = [];

		for (let i = 0; i < paths.length; i++) {
			if (paths[i].length < sequence.length) continue;

			paths[i].push(
				paths[i][sequence.length - 1].dvused + 
				(includecapture ? this.calcEjectiondV(
					sequence[sequence.length - 1],
					paths[i][sequence.length - 1].prevtransfer.incomingmag,
					truefinalalt
				) : 0)
			);

			let duration = paths[i][sequence.length - 1].prevtransfer.t2 - paths[i][1].prevtransfer.t1;

			paths[i].push(duration);

			if (duration > maxduration) continue;

			trajectories.push(paths[i]);
		}

		trajectories.sort(function(a, b) {
			return a[sequence.length] - b[sequence.length];
		});

		if (trajectories.length === 0) return [];

		console.log(trajectories);
		console.log(trajectories[0]);

		// Local optimization to optimize best trajectory
		// For now, I'll use a simple hill climb but in the future I might use a different algorithm here.
		let besttrajectory = [];

		for (let i = 1; i < sequence.length; i++) {
			let t1 = trajectories[0][i].prevtransfer.t1;
			let t2 = trajectories[0][i].prevtransfer.t2;

			besttrajectory.push(t1);
			besttrajectory.push(trajectories[0][i].prevtransfer.outgoing.x);
			besttrajectory.push(trajectories[0][i].prevtransfer.outgoing.y);
			besttrajectory.push(trajectories[0][i].prevtransfer.outgoing.z);

			if (trajectories[0][i].prevtransfer.hasOwnProperty("tdsm")) {
				besttrajectory.push((trajectories[0][i].prevtransfer.tdsm - t1) / (t2 - t1));
			} else {
				besttrajectory.push(0.5);
			}
		}

		besttrajectory.push(trajectories[0][sequence.length - 1].prevtransfer.t2);
		let bestfitness = this.fitness(besttrajectory, sequence, trueinitalt, truefinalalt, includecapture);

		console.log(besttrajectory);
		console.log(bestfitness);

		for (let i = 0; i < 10000; i++) {
			let newtrajectory = structuredClone(besttrajectory);

			let curr = 0;
			for (let j = 0; j < sequence.length; j++) {
				newtrajectory[5 * j] += normal() * 100;

				if (newtrajectory[5 * j] < curr) newtrajectory[5 * j] = curr + 1;

				curr = newtrajectory[j];
			}

			for (let j = 0; j < sequence.length - 1; j++) {
				newtrajectory[5 * j + 1] += normal() * 10;
				newtrajectory[5 * j + 2] += normal() * 10;
				newtrajectory[5 * j + 3] += normal() * 10;
				newtrajectory[5 * j + 4] += normal() * 0.01;

				if (newtrajectory[5 * j + 4] < 0.1) {
					newtrajectory[5 * j + 4] = 0.1;
				}
				if (newtrajectory[5 * j + 4] > 0.9) {
					newtrajectory[5 * j + 4] = 0.9;
				}
			}

			let newfitness = this.fitness(newtrajectory, sequence, trueinitalt, truefinalalt, includecapture);

			if (newfitness < bestfitness) {
				bestfitness = newfitness;
				besttrajectory = structuredClone(newtrajectory);

				console.log(newtrajectory, newfitness);
			}
		}

		return [besttrajectory, bestfitness];
	}
};

let mgafinder;

onmessage = (e) => {
	let copiedsystem = e.data.system;
	for (let i = 0; i < copiedsystem.bodies.length; i++) {
		copiedsystem.bodies[i].orbit = Object.setPrototypeOf(copiedsystem.bodies[i].orbit, Orbit.prototype);
	}

	mgafinder = new MGAFinder(Object.setPrototypeOf(copiedsystem, System.prototype));

	let params = e.data.params;

	console.log(params);

	let res = mgafinder.planMGATrajectory(...params);

	if (res.length === 0) {
		postMessage({
			status: "failure"
		});
	} else {
		res.push(params);

		postMessage({
			status: "success",
			result: res
		});
	}
}