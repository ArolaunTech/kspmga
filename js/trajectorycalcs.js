import { Vector3 } from './vector.js';
import { brentMinimize } from './optimize.js';

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

export function calcEjectionDetailsInclined(body, vout, circalt, system) {
	let vinf = vout.norm;
	let hvinf = Math.hypot(vout.x, vout.y);
	let vvinf = vout.z;
	let rotvout = new Vector3(hvinf, 0, vvinf);

	let stdgp = system.bodies[system.bodymap.get(body)].gravparameter;
	let soi = system.bodies[system.bodymap.get(body)].soi;
	let radius = system.bodies[system.bodymap.get(body)].radius;
	let minalt = radius + system.bodies[system.bodymap.get(body)].atmodepth + 10000;

	let realcircalt = circalt + radius;

	let energy = 0.5 * vinf * vinf;
	let sma = -stdgp / 2 / energy;

	let circvel = Math.sqrt(stdgp / realcircalt);
	let ejectionvel = Math.sqrt(2 * (energy + stdgp / realcircalt));

	let xdir = Vector3.normalize(rotvout);

	let x = brentMinimize(
		(x)=>{
			let R = new Vector3(realcircalt * Math.cos(x), realcircalt * Math.sin(x), 0);
			let ydir = Vector3.normalize(Vector3.cross(Vector3.cross(xdir, R), xdir));
			let b = 0.5 * (Vector3.dot(ydir, R) + Math.sqrt(Math.pow(Vector3.dot(ydir, R), 2) + 4 * sma * (Vector3.dot(xdir, R) - realcircalt)));
			let h = b * vinf;
			let p = h * h / stdgp;
			let hvec = Vector3.mult(Vector3.normalize(Vector3.cross(R, rotvout)), h);
			let evec = Vector3.sub(Vector3.mult(Vector3.cross(rotvout, hvec), 1 / stdgp), Vector3.normalize(rotvout));
			let v0 = Vector3.mult(Vector3.cross(hvec, Vector3.add(evec, Vector3.normalize(R))), 1 / p);
			let vprev = new Vector3(-circvel * Math.sin(x), circvel * Math.cos(x), 0);
			let dv = Vector3.sub(v0, vprev).norm;

			let e = Math.sqrt(1 - p / sma);
			let periapsis = sma * (1 - e);

			let rvdot = Vector3.dot(v0, R);

			if (rvdot > 0 || periapsis > minalt) return dv;
			return dv + 1e6 * (minalt - periapsis); // fly safe!
		},
		Math.PI,
		Math.PI * 2,
		1e-8,
		50
	)[0];

	let R = new Vector3(realcircalt * Math.cos(x), realcircalt * Math.sin(x), 0);
	let ydir = Vector3.normalize(Vector3.cross(Vector3.cross(xdir, R), xdir));
	let b = 0.5 * (Vector3.dot(ydir, R) + Math.sqrt(Math.pow(Vector3.dot(ydir, R), 2) + 4 * sma * (Vector3.dot(xdir, R) - realcircalt)));
	let h = b * vinf;
	let p = h * h / stdgp;
	let hvec = Vector3.mult(Vector3.normalize(Vector3.cross(R, rotvout)), h);
	let evec = Vector3.sub(Vector3.mult(Vector3.cross(rotvout, hvec), 1 / stdgp), Vector3.normalize(rotvout));
	let v0 = Vector3.mult(Vector3.cross(hvec, Vector3.add(evec, Vector3.normalize(R))), 1 / p);
	let vprev = new Vector3(-circvel * Math.sin(x), circvel * Math.cos(x), 0);
	let maneuver = Vector3.sub(v0, vprev);

	let prograde = Vector3.dot(maneuver, new Vector3(-Math.sin(x), Math.cos(x), 0));
	let radial = Vector3.dot(maneuver, new Vector3(Math.cos(x), Math.sin(x), 0));
	let normal = Vector3.dot(maneuver, new Vector3(0, 0, 1));

	let rvdot = Vector3.dot(v0, R);
	let e = Math.sqrt(1 - p / sma);
	let c0 = e * b * Math.sqrt(stdgp / p);
	let F0 = Math.asinh(rvdot / c0);
	let M0 = e * rvdot / c0 - F0;

	let F1 = Math.acosh((1 - soi / sma) / e);
	let M1 = e * Math.sinh(F1) - F1;

	let dt = sma * Math.sqrt(-sma / stdgp) * (M0 - M1);

	return {
		prograde: prograde,
		radial: radial,
		normal: normal,
		dv: Math.hypot(prograde, radial, normal),
		dt: dt
	};
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