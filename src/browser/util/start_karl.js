//-----------------------------------------------------------------------------------------------------
// Initialize the global variable karl. All the physics routines are functions that are members of this
// variable.
//-----------------------------------------------------------------------------------------------------
var karl = {};
karl.load = function(){}; /* Used in the rhino implementation. In the browser implementation, we load
                             scripts using the html <script> tag. */

//-----------------------------------------------------------------------------------------------------
// Provide an implementation of the rhino print() function. If you want this output to appear on
// the web page, make a scrollable object called "console", e.g.:
//   <div style="height:300px;width:800px;overflow:auto" id="console"></div>
// As with rhino's print(), it can be called with any number of arguments, and those get converted to strings.
//-----------------------------------------------------------------------------------------------------
var verbosity = 1; // Debugging output gets turned on if this is 2 or 3.
function print() {
  var args = [];
  for (var i=0; i<arguments.length; i++) {
    args.push(arguments[i]);
    if (i<arguments.length-1) {
      args.push(' '); // emulate the behavior of my python implementation
    }
  }
  var m = io_util.strcat(args);
  var output_area = document.getElementById("console");
  var t = output_area.innerHTML;
  if (t!="") {t=t+"<br>"}
  output_area.innerHTML = t + m;
  // Make it scroll down automatically. See https://stackoverflow.com/a/21067431/1142217
  var isScrolledToBottom = output_area.scrollHeight - output_area.clientHeight <= output_area.scrollTop + 1;
  if (!isScrolledToBottom) {output_area.scrollTop = output_area.scrollHeight - output_area.clientHeight;}
  
}
