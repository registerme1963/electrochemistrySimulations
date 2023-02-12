Array.prototype.Transpose = function() { var a = this; return Object.keys(a[0]).map(function (c) { return a.map(function (r) { return r[c]; }); });	}
Array.prototype.Mean = function() { var x = this, sum = x.reduce(function(a, b) { return a + b; }); return sum / x.length; }

function roundToTwo(num) {
    return +(Math.round(num + "e+2")  + "e-2");
}

var DBG = true;

function mylog(str) {
	if (DBG) { console.log(str); }
	else {  }
}

function randomKey() {
    return Math.random().toString(36).substring(2, 12);
}

/************************
 * AngularJS code
 */

var app = angular.module('app', ['ngAria', 'ngMessages', 'ngAnimate', 'ngMaterial', 'ngSanitize', 'chart.js']);

app.controller('ECSimCtrl', function($scope, $timeout, $mdDialog) {

	// global variables:
	$scope.selected_tab = 0; // start at the DATA tab
	
	// simulation system and results:
	$scope.simsys = {};
    $scope.simoutput = "";
    $scope.simoutputclipboard = "";
    $scope.expinputclipboard = "";
    
    
    // species:
    $scope.simsys.species = [];
    $scope.simsys.species.push( {id: 1, name: "ox1", diff: 1.0e-9, conc: 1.0} );
    $scope.simsys.species.push( {id: 2, name: "ox2", diff: 1.0e-9, conc: 0.0} );
    $scope.simsys.species.push( {id: 3, name: "red1", diff: 1.0e-9, conc: 0.0} );
    $scope.simsys.species.push( {id: 4, name: "red2", diff: 1.0e-9, conc: 0.0} );
    $scope.simsys.species.push( {id: 5, name: "prod", diff: 1.0e-9, conc: 0.0} );
    // AddSpecies adds a species to $scope.simsys.species array with a random name
    $scope.AddSpecies = function() {
        var maxSpeciesId = Math.max.apply(Math, $scope.simsys.species.map(function(s) { return s.id; }));
        $scope.simsys.species.push( {id: maxSpeciesId+1, name: randomKey(), diff: 1.0e-9, conc: 0.0} );
    };
    // RemoveSpecies(id) removes species from $scope.simsys.species array (if it is not in use)
    $scope.RemoveSpecies = function(id) {
        var index = $scope.simsys.species.map(function(s) { return s.id; }).indexOf(id);
        if (!$scope.speciesIsInUse(id) && index >= 0) $scope.simsys.species.splice(index, 1);
    };
    // speciesIsInUse(id) returns true if species is in use in redox or reaction
    $scope.speciesIsInUse = function(id) {
        var isSpeciesInUse = false;
        angular.forEach($scope.simsys.redox, function(r) { isSpeciesInUse |= (r.ox == id || r.red == id); });
        angular.forEach($scope.simsys.rxn,   function(c) { isSpeciesInUse |= (c.reactant1 == id || c.reactant2 == id || c.product1 == id || c.product2 == id); });
        return isSpeciesInUse;
    };
    // speciesIdToName(id) returns the name of the species if it exists, otherwise empty string
    $scope.speciesIdToName = function(id) {
        var name = "";
        angular.forEach($scope.simsys.species, function(s) { if (s.id == id) name = s.name; });
        return name;
    };
    
    // redox:
    $scope.simsys.redox = [];
    $scope.simsys.redox.push( {ox: 1, red: 3, ne: 1, pot: -0.5, ke: 1.0, alpha: 0.5, enabled: 1} );
    $scope.simsys.redox.push( {ox: 2, red: 4, ne: 1, pot: 0.5, ke: 1.0e-5, alpha: 0.5, enabled: 1} );
    $scope.AddRedox = function() { $scope.simsys.redox.push( {ox: "", red: "", ne: 1, pot: -0.5, ke: 1.0, alpha: 0.5, enabled: 0} ); };
    $scope.RemoveRedox = function(idx) { $scope.simsys.redox.splice(idx, 1); };
    // reaction:
    $scope.simsys.rxn = [];
    $scope.simsys.rxn.push( {reactant1: 3, reactant2: 0, product1: 4, product2: 0, kf: 10.0, kb: 0.0, enabled: 1} );
    $scope.simsys.rxn.push( {reactant1: 3, reactant2: 2, product1: 5, product2: 0, kf: 10.0, kb: 0.0, enabled: 0} );
    $scope.AddReaction = function() { $scope.simsys.rxn.push( {reactant1: 0, reactant2: 0, product1: 0, product2: 0, kf: 10.0, kb: 0.0, enabled: 0} ); };
    $scope.RemoveReaction = function(idx) { $scope.simsys.rxn.splice(idx, 1); };
    
    
    // electrode:
    $scope.simsys.electrodetypes = [];
    $scope.simsys.electrodetypes.push( {type: 0, typename: "Disk", geom1name: "Radius", geom2name: ""} );
    $scope.simsys.electrodetypes.push( {type: 1, typename: "Square", geom1name: "Width", geom2name: ""} );
    $scope.simsys.electrodetypes.push( {type: 2, typename: "Rectangle", geom1name: "Width", geom2name: "Length"} );
    $scope.simsys.electrodetypes.push( {type: 3, typename: "Cylinder", geom1name: "Radius", geom2name: "Length"} );
    $scope.simsys.electrodetypes.push( {type: 4, typename: "Sphere", geom1name: "Radius", geom2name: ""} );
    $scope.simsys.electrodetypes.push( {type: 5, typename: "Hemisphere", geom1name: "Radius", geom2name: ""} );
    $scope.simsys.electrode = {type: 0, geom1: 0.001, geom2: 0.001};
    
    // environment:
    $scope.simsys.temp = 293.15;
    
    // experiment:
    $scope.experiment = {};
    $scope.experiment.conditioningPotential = 0.0;
    $scope.experiment.conditioningTime = 0.0;
    $scope.experiment.equilibrationTime = 0.0;
    
    $scope.experiment.initialPotential = 0.0;
    $scope.experiment.vertexPotentials = [ {pot: -0.8, key: randomKey()}, {pot: 0.9, key: randomKey()} ]; // randomkey necessary to avoid "dupes" error on ng-repeat
    $scope.AddVertex = function() { $scope.experiment.vertexPotentials.push( {pot: 0.0, key: randomKey()} ); };
    $scope.RemoveVertex = function(idx) { $scope.experiment.vertexPotentials.splice(idx, 1); };
    $scope.experiment.finalPotential = 0.0;
    
    $scope.experiment.vertexDelay = 0.0;
    $scope.experiment.numCycles = 1;
    $scope.experiment.scanRate = 1.0;
    
    // simulation settings:
    $scope.simsettings = { Nderivative: 6, Ncurrent: 5, deltatheta: 0.2, gammae: 1.67, F0: 2.2, F1: 6.2, lograte0: 3.0, lograte1: 7.0 };
    
    // results:
    $scope.color_palette = palette('tol', 10).map(c => '#' + c);
    $scope.chartsim = {
		series: [],
		data: [ [], [] ], // [0] is simulated data, [1] is experimental data from clipboard
		colors: [ $scope.color_palette[0] ],
		options: {
			animation: false,
			scales: {
				xAxes: [{ type: 'linear', position: 'bottom', scaleLabel: { display: true, labelString: 'Potential [V]' } }],
				yAxes: [{ type: 'linear', position: 'left',   scaleLabel: { display: true, labelString: 'Current [A]' } }]
			},
			legend: { display: false, position: 'top' }
		},
		datasetOverride: []
	};
    
    $scope.RunSimulation = function() {
        
        console.log( $scope.simsys );
        
        // clear simulation:
        Module.ccall("clearSim", "number", [], []);
    
        // set species:
        angular.forEach($scope.simsys.species, function(s) {
            Module.ccall("addSpecies", "number",
                ["string", "number", "number"],
                [s.name, s.conc, s.diff]);
        });
        
        // set redox:
        angular.forEach($scope.simsys.redox, function(r) {
            if (r.ox > 0 && r.red > 0 && r.enabled) {
                
                Module.ccall("addRedox", "number",
                    ["string", "string", "number", "number", "number", "number", "number", "number"],
                    [$scope.speciesIdToName(r.ox), $scope.speciesIdToName(r.red),
                             r.ne, r.pot, r.ke, r.alpha, 0, r.enabled]);
            }
        });
        
        // set reactions:
        angular.forEach($scope.simsys.rxn, function(c) {
            if (c.reactant1 > 0 && c.product1 > 0 && c.enabled)
                Module.ccall("addReaction", "number",
                    ["string", "string", "string", "string", "number", "number", "number"],
                    [$scope.speciesIdToName(c.reactant1), $scope.speciesIdToName(c.reactant2),
                             $scope.speciesIdToName(c.product1), $scope.speciesIdToName(c.product2),
                             c.kf, c.kb, c.enabled]);
        });
        
        // set electrode:
        Module.ccall("setElectrode", "number",
            ["number", "number", "number"],
            [$scope.simsys.electrode.type, $scope.simsys.electrode.geom1, $scope.simsys.electrode.geom2]);
        
        // set environment:
        Module.ccall("setEnvironment", "number", ["number"], [$scope.simsys.temp]);
        
        // set experiment:
        Module.ccall("setCVExperiment", "number",
            ["number", "number", "number", "number", "number", "number", "number", "number"],
            [$scope.experiment.conditioningTime, $scope.experiment.conditioningPotential, $scope.experiment.equilibrationTime, $scope.experiment.initialPotential,
                $scope.experiment.finalPotential, $scope.experiment.vertexDelay, $scope.experiment.scanRate, $scope.experiment.numCycles]);
        angular.forEach($scope.experiment.vertexPotentials, function(v) { Module.ccall("addCVVertex", "number", ["number"], [v.pot]); });
        Module.ccall("setCVVertices", "number", [], []);
        
        // set simulation settings:
        Module.ccall("setGridSizing", "number",
            ["number", "number", "number", "number", "number"],
            [$scope.simsettings.gammae, $scope.simsettings.F0, $scope.simsettings.F1, $scope.simsettings.lograte0, $scope.simsettings.lograte1]);
        Module.ccall("setPotentialSizing", "number", ["number"], [$scope.simsettings.deltatheta]);
        Module.ccall("setDifferentialOrders", "number", ["number", "number"], [$scope.simsettings.Ncurrent, $scope.simsettings.Nderivative]);
        
        // run simulation:
        Module.simOutput = "";
        var sz = Module.ccall('dosim', 'number', [], []); // run simulation; returns number of points in simulation
        var ptr_current = Module.ccall('getcurrent', 'number', [], []); // get pointer to current array on heap
        var ptr_potential = Module.ccall('getpotential', 'number', [], []); // get pointer to potential array on heap
        var points = []; // array of [potential, current] pairs
        
        $scope.chartsim.data[0] = [];
        $scope.simoutputclipboard = "";
        for (var i = 0; i < sz; i++) {
            var potval =  Module.HEAPF64[ptr_potential/Float64Array.BYTES_PER_ELEMENT+i];
            var currval = Module.HEAPF64[ptr_current/Float64Array.BYTES_PER_ELEMENT+i];
            $scope.chartsim.data[0].push({x: potval, y: currval});
            $scope.simoutputclipboard += potval.toString() + "\t" + currval.toString() + "\r\n";
        }
        
        console.log( Module.simOutput );
        $scope.simoutput = Module.simOutput;
    };
    
    $scope.DataPasted = function($event) {
        if(typeof $event.originalEvent == "undefined") {
            pastedData = event.clipboardData.getData('text/plain');
        } else if(typeof $event.originalEvent.clipboardData !== "undefined") {
            pastedData = $event.originalEvent.clipboardData.getData('text/plain');
        } else { // To support browsers without clipboard API (IE and older browsers)
            $timeout(function(){
                pastedData = angular.element($event.currentTarget).val();
            });
        }
        
        console.log(pastedData);

        var rowStrings = pastedData.split(/[\r\n]+/);
        var numRows = rowStrings.length;
        var trimmedRowString = "";
        var colStrings = "";
        
        $scope.chartsim.data[1] = [];
        for (var r = 0; r < numRows; r++) {
            trimmedRowString = rowStrings[r].trim();
            if (trimmedRowString) { // check if not empty or undefined
                colStrings = trimmedRowString.split(/,?[\t\s]+/);
                if (colStrings.length == 2) {
                    $scope.chartsim.data[1].push( {x: parseFloat(colStrings[0]), y: parseFloat(colStrings[1])} );
                }
            }
        }
        
        $scope.expinputclipboard = "";
    }
    
    $scope.CopyDataToClipboard = function() {
        // stolen from: https://stackoverflow.com/questions/46041831/copy-to-clipboard-with-break-line
        const myFluffyTextarea = document.createElement('textarea');    
        myFluffyTextarea.innerHTML = $scope.simoutputclipboard;
        const parentElement = document.getElementById('simresults');
        parentElement.appendChild(myFluffyTextarea);
        myFluffyTextarea.select();
        document.execCommand('copy');
        parentElement.removeChild(myFluffyTextarea);
    };
    
	function DialogController($scope, $mdDialog) {
		$scope.hide = function() { $mdDialog.hide(); };
		$scope.cancel = function() { $mdDialog.cancel(); };
		$scope.answer = function(answer) { $mdDialog.hide(answer); };
	}
	$scope.UserMessage = function(type, message) {
		if (type == 'alert') {
			$mdDialog.show( $mdDialog.alert()
				.clickOutsideToClose(true)
				.title(message[0])
				.htmlContent(message[1])
				.ariaLabel('Alert dialog')
				.ok('Got it!')
			);
		} else if (type == 'info') {
			console.log(message);
			$mdDialog.show({
				controller: DialogController,
				templateUrl: message,
				parent: angular.element(document.body),
				clickOutsideToClose:true,
				fullscreen: false
			});
		}
	};
});
