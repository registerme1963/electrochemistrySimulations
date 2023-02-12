

/************************
 * WebAssembly code
 */

function htmlDecode(input)
{
  var doc = new DOMParser().parseFromString(input, "text/html");
  return doc.documentElement.textContent;
}

var Module = {
    simOutput: "",
    preRun: [],
    postRun: [],
    print: (function() {
        return function(text) {
        if (arguments.length > 1) text = Array.prototype.slice.call(arguments).join(' ');
        //outputElement.append(text + "\n");
        Module.simOutput += htmlDecode(text) + "\n";
        };
    })(),
    printErr: function(text) {
        if (arguments.length > 1) text = Array.prototype.slice.call(arguments).join(' ');
        console.error(text);
    },
    canvas: (function() {
        return document.getElementById("wasmcanvas");
    })(),
    setStatus: function(text) { console.log(text); },
    totalDependencies: 0,
    monitorRunDependencies: function(left) {
        this.totalDependencies = Math.max(this.totalDependencies, left);
        Module.setStatus(left ? 'Preparing... (' + (this.totalDependencies-left) + '/' + this.totalDependencies + ')' : 'All downloads complete.');
    }
};
Module.setStatus('Downloading...');
window.onerror = function(event) {
    // TODO: do not warn on ok events like simulating an infinite loop or exitStatus
    Module.setStatus('Exception thrown, see JavaScript console');
    Module.setStatus = function(text) {
        if (text) Module.printErr('[post-exception status] ' + text);
    };
};

