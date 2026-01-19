import { Vector3 } from './vector.js';

export function calcEjectiondV(body, vinf, circalt, system) {
	let stdgp = system.bodies[system.bodymap.get(body)].gravparameter;
	let soi = system.bodies[system.bodymap.get(body)].soi;
	let radius = system.bodies[system.bodymap.get(body)].radius;

	let energy = 0.5 * vinf * vinf - stdgp / soi;
	let sma = -stdgp / 2 / energy;

	let ejectionvel = Math.sqrt(stdgp * (2 / (circalt + radius) - 1 / sma));
	let circvel = Math.sqrt(stdgp / (circalt + radius));

	return ejectionvel - circvel;
}

export function calcPeriapsis(body, vin, vout, system) {
	let stdgp = system.bodies[system.bodymap.get(body)].gravparameter;
	let radius = system.bodies[system.bodymap.get(body)].radius;
	let minalt = system.bodies[system.bodymap.get(body)].atmodepth + 10000;

	//console.log(stdgp, radius, minalt, vin, vout, body);

	let cosdeflection = Vector3.dot(vin, vout) / vin.norm / vout.norm;
	let rp = (1 / Math.sqrt((1 - cosdeflection) / 2) - 1) * stdgp / vin.norm / vin.norm;

	//console.log(cosdeflection, rp);

	rp -= radius;

	if (rp < minalt) rp = minalt;

	//console.log(rp);

	return rp;
}

export function calcCosDeflection(body, vinf, system) {
	let stdgp = system.bodies[system.bodymap.get(body)].gravparameter;
	let minalt = system.bodies[system.bodymap.get(body)].radius + system.bodies[system.bodymap.get(body)].atmodepth + 10000;

	return 1 - 2 / Math.pow(1 + minalt * vinf * vinf / stdgp, 2);
}

export function violation(body, vin, vout, system) {
	let innorm = vin.norm;
	let outnorm = vout.norm;

	let cosmaxdeflection = calcCosDeflection(body, innorm, system);
	let cosdeflection = Vector3.dot(vin, vout) / innorm / outnorm;

	if (cosdeflection > cosmaxdeflection) {
		return Math.abs(outnorm - innorm);
	}

	let sindir = Vector3.normalize(Vector3.cross(Vector3.cross(vin, vout), vin));
	let closest = Vector3.add(Vector3.mult(vin, cosmaxdeflection), Vector3.mult(sindir, innorm * Math.sqrt(1 - cosmaxdeflection * cosmaxdeflection)));

	return Vector3.sub(closest, vout).norm;
}