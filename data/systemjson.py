import os
import pyparsing as pp
import json
import math

#========== Inputs ==========
folderpath = "configs/stock"
outputpath = "systems/stock.json"
systemname = "stock"

#========== Consts ==========
lbracket = pp.Literal("{").suppress()
rbracket = pp.Literal("}").suppress()
equals = pp.Literal("=").suppress()

whitespace = pp.White()[...].suppress()

label = pp.Word(pp.printables, excludeChars="{}")
equation = (label + equals + whitespace + pp.SkipTo(pp.LineEnd() | "}")).set_name('equation')

node = pp.Forward()

atom = pp.Group(equation) | node
node <<= pp.Group(label + lbracket + atom[...] + rbracket)

node.set_name('node')

config = pp.OneOrMore(atom).set_name('config')
config.create_diagram("out.html", show_hidden=True)

configs = os.listdir(folderpath)

res = {"name": "stock", "bodies": []}

rootmass = -1
for configpath in configs:
	if configpath[-3:] != "cfg":
		continue

	configstring = ""
	with open(f"{folderpath}/{configpath}", "r") as f:
		for line in f:
			newline = line
			if "//" in newline:
				newline = newline[:newline.index("//")] + "\n"
			if newline == "\n":
				continue

			configstring += newline

	parsed = config.parse_string(configstring)[0]
	bodynode = parsed[1]
	atmonode = []
	orbitnode = []
	
	for eq in bodynode:
		if eq[0] == 'Properties':
			propertynode = eq
		if eq[0] == "Atmosphere":
			atmonode = eq
		if eq[0] == "Orbit":
			orbitnode = eq
		if eq[0] == 'name':
			bodyname = eq[1]

	for eq in propertynode:
		if eq[0] == 'radius':
			radius = float(eq[1])

	soi = -1
	for eq in propertynode:
		if eq[0] == 'mass':
			gravparameter = 6.68408e-11 * float(eq[1])
		if eq[0] == 'gravParameter':
			gravparameter = float(eq[1])
		if eq[0] == 'geeASL':
			gravparameter = float(eq[1]) * radius * radius * 9.80665
		if eq[0] == 'sphereOfInfluence':
			soi = float(eq[1])

	atmosphere = False
	atmoalt = -1
	if atmonode != []:
		for eq in atmonode:
			if eq[0] == 'enabled':
				atmosphere = (eq[1].lower() == 'true')
			if eq[0] in ['atmosphereDepth', 'altitude', 'maxAltitude']:
				atmoalt = float(eq[1])

	root = False
	sma = -1
	eccentricity = -1
	inclination = -1
	lan = -1
	argp = -1
	initanomaly = -1
	epoch = 0
	refbody = None
	color = [1, 1, 0.75, 1]
	if orbitnode == []:
		root = True
		rootmass = gravparameter
	else:
		for eq in orbitnode:
			if eq[0] == "semiMajorAxis":
				sma = float(eq[1])
			if eq[0] == "eccentricity":
				eccentricity = float(eq[1])
			if eq[0] == 'inclination':
				inclination = float(eq[1]) * math.pi / 180
			if eq[0] == "referenceBody":
				refbody = eq[1]
			if eq[0] == "longitudeOfAscendingNode":
				lan = float(eq[1]) * math.pi / 180
			if eq[0] == "argumentOfPeriapsis":
				argp = float(eq[1]) * math.pi / 180
			if eq[0] == "meanAnomalyAtEpoch":
				initanomaly = float(eq[1])
			if eq[0] == "meanAnomalyAtEpochD":
				initanomaly = float(eq[1]) * math.pi / 180
			if eq[0] == "epoch":
				epoch = float(eq[1])
			if eq[0] == "color":
				# I only support Fractional RGBA lists right now
				rgba = eq[1].split(",")
				color = list(map(float, rgba))

	res["bodies"].append({
		"name": bodyname, 
		"radius": radius, 
		"gravparameter": gravparameter,
		"atmospheric": atmosphere,
		"atmodepth": atmoalt,
		"root": root,
		"soi": soi,
		"color": color,
		"orbit": {
			"sma": sma,
			"eccentricity": eccentricity,
			"inclination": inclination,
			"parent": refbody,
			"lan": lan,
			"argp": argp,
			"meananomalyatepoch": initanomaly,
			"epoch": epoch
		}
	})

for i in range(len(res["bodies"])):
	if res["bodies"][i]["soi"] == -1 and not res["bodies"][i]["root"]:
		res["bodies"][i]["soi"] = res["bodies"][i]["orbit"]["sma"] * (res["bodies"][i]["gravparameter"] / rootmass) ** 0.4
	if res["bodies"][i]["root"]:
		res["bodies"][i]["soi"] == 1e1000

with open(outputpath, "w") as f:
	f.write(json.dumps(res, indent=4).replace("Infinity", "1e1000"))