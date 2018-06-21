PJ and translate_maxima
=======================

## Purpose

These are two small utilities, packaged with Karl, that help with translation to javascript.

## PJ

The simplest, stupidest python-to-javascript converter that has any chance of working.
Limitations: 

* fine-tuned for my application, not for general-purpose use
* doesn't handle any OOP
* needs explicit \ characters for continuation lines in all cases
* is so simple that certain stuff will always have to be translated by hand

Hand-translation is handled using active comments that look like this:

    x = [0 for i in range(2)] \
    #js var x=[0,0];

(The active comment does not have to be on a continuation line if you'd rather have it on the same line.)

For js code that doesn't have any counterpart in the python code:

    #js foo(bar);

In assignment statements, if certain keywords such as sin, log, **, ... occur in the
rhs of the source code, run the rhs through translate_maxima.

The script should be invoked with a command-line argument that specifies the name of the module, to
be prepended to all function definitions. Following python's convention, this would normally be the
same as the name of the file. As a convenience, this command-line argument can be the full filename,
of which only the stem of the filename will be used. So a typical invocation of the script would be:

    pj.rb foo.py <foo.py >foo.js

### Handling of modules and loading

The translator looks at the top of a source-code file for a #! line in order to tell
whether it's a module or the main program. At runtime, the generated code also uses
the IS_BROWSER macro to determine whether it's running in a browser or in an environment
like rhino from the command line.

The main program automatically has a header generated, which initializes a global
called karl, and a hash karl.modules_loaded. Import statements in python translate
to code which, in a non-browser environment, causes a module to be loaded (through
lib/loader.js) using
rhino's load() method, and the module is marked in karl.modules_loaded so that
it won't inadvertently be loaded twice. 

A module also automatically has a header generated. For a module named foo, this
header causes a global variable called foo to be created, and if the module defines
a function named bar, then a member foo.bar is created as a reference to that function.
Thus python syntax like foo.bar(x,y) is also valid in javascript.

The module-loading stuff is not used for the modules test.js, lib/math.js, and lib/lambertw.js.
It should be, but I coded it after I set up those modules.

### Not yet implemented, may not need

To control the behavior of the translator, can do active comments like this:

    #js:{"foo":2}

The stuff after the : is a JSON hash. In this example, the foo parameter is set to 2.

## translate_maxima

This translates a mathematical expression written in the computer algebra system Maxima
into javascript. Reads stdin, writes stdout. Example:

    echo "sin(x)+x^y+whacko(p,q)/2" | ./translate_maxima.py
    ((((((1.0)/(2.0)))*(whacko((p),(q)))))+(math.sin((x)))+(math.pow((x),(y))))

The output language doesn't have to be javascript, could be pretty much any language with
infix notation. The language-specific stuff is controlled through some global data
structures at the top of the script.

If there's an error, then the output line consists of a caret followed by the error message.

