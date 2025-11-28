function secsToKerbalDuration(secs) {
	let seconds = secs % 60;
	let minutes = Math.floor(secs / 60) % 60;
	let hours = Math.floor(secs / 3600) % 6;
	let days = Math.floor(secs / 21600) % 426;
	let years = Math.floor(secs / 9201600);

	return [years, days, hours, minutes, seconds];
}

function secsToKerbalTime(secs) {
	let duration = secsToKerbalDuration(secs);

	duration[0]++;
	duration[1]++;

	return duration;
}

function secsToKerbalTimeString(secs) {
	let time = secsToKerbalTime(secs);

	let stringminutes = String(time[3]);
	let stringseconds = String(time[4]);

	if (stringminutes.length === 1) {
		stringminutes = '0' + stringminutes;
	}
	if (stringseconds.length === 1) {
		stringseconds = '0' + stringseconds;
	}

	return `Year ${time[0]}, Day ${time[1]} ${time[2]}:${stringminutes}:${stringseconds}`;
}

function kerbalTimeToSecs(time) {
	let duration = time;

	duration[0]--;
	duration[1]--;

	return kerbalDurationToSecs(duration);
}

function kerbalDurationToSecs(duration) {
	return 9201600 * duration[0] + 21600 * duration[1] + 3600 * duration[2] + 60 * duration[3] + duration[4];
}