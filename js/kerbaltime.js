var SECONDSPERDAY = 21600;
var SECONDSPERYEAR = 9201600;

export function settimesystem(secsperday, secsperyear) {
	SECONDSPERDAY = secsperday;
	SECONDSPERYEAR = secsperyear;
}

export function gettimesystem() {
	return [SECONDSPERDAY, SECONDSPERYEAR];
}

export function secsToKerbalDuration(secs) {
	let seconds = secs % 60;
	let minutes = Math.floor(secs / 60) % 60;
	let hours = Math.floor(secs / 3600) % 6;
	let days = Math.floor(secs / SECONDSPERDAY) % (SECONDSPERYEAR / SECONDSPERDAY);
	let years = Math.floor(secs / SECONDSPERYEAR);

	return [years, days, hours, minutes, seconds];
}

export function secsToKerbalTime(secs) {
	let duration = secsToKerbalDuration(secs);

	duration[0]++;
	duration[1]++;

	return duration;
}

export function secsToKerbalTimeString(secs) {
	let time = secsToKerbalTime(secs);

	let stringminutes = String(time[3]);
	let stringseconds = String(Math.round(time[4]));

	if (stringminutes.length === 1) {
		stringminutes = '0' + stringminutes;
	}
	if (stringseconds.length === 1) {
		stringseconds = '0' + stringseconds;
	}

	return `Year ${time[0]}, Day ${time[1]} ${time[2]}:${stringminutes}:${stringseconds}`;
}

export function kerbalTimeToSecs(time) {
	let duration = time;

	duration[0]--;
	duration[1]--;

	return kerbalDurationToSecs(duration);
}

export function kerbalDurationToSecs(duration) {
	return SECONDSPERYEAR * duration[0] + SECONDSPERDAY * duration[1] + 3600 * duration[2] + 60 * duration[3] + duration[4];
}