PJ and translate_maxima
=======================

## Purpose

These are two small utilities, packaged with Karl, that help with translation to javascript.

## translate_maxima

This translates a mathematical expression written in the computer algebra system Maxima
into javascript. Reads stdin, writes stdout. Example:

    echo "sin(x)+x^y+whacko(p,q)/2" | ./translate_maxima.py
    ((((((1.0)/(2.0)))*(whacko((p),(q)))))+(math.sin((x)))+(math.pow((x),(y))))

The output language doesn't have to be javascript, could be pretty much any language with
infix notation. The language-specific stuff is controlled through some global data
structures at the top of the script.

If there's an error, then the output line consists of a caret followed by the error message.

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

### Not yet implemented, may not need

To control the behavior of the translator, can do active comments like this:

    #js:{"foo":2}

The stuff after the : is a JSON hash. In this example, the foo parameter is set to 2.

