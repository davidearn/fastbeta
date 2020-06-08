![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") in R package documentation with mathjaxr
================
Mikael Jagan
07/06/2020

This document contains evolving notes on the
[mathjaxr](https://github.com/vnijs/MathJaxR) package, developed by
Wolfgang Viechtbauer in response to mailing list threads
[1](https://stat.ethz.ch/pipermail/r-devel/2020-April/079395.html) and
[2](https://stat.ethz.ch/pipermail/r-package-devel/2020q2/005420.html).
The package automates integration of
[MathJax](https://en.wikipedia.org/wiki/MathJax) in R package
documentation, making it much easier to get typeset math in HTML help
pages.

## Installation

The current development version of mathjaxr (1.0-0) is much more useful
than current CRAN version (0.8-3). This will change with the next
update, so check
[here](https://cran.r-project.org/web/packages/mathjaxr/index.html)
before installing.

``` r
# From CRAN
install.packages("mathjaxr")

# From GitHub (recommended until release of 1.0-0)
remotes::install_github("wviechtb/mathjaxr")
```

## Configuring `DESCRIPTION` and `NAMESPACE`

The mathjaxr package defines macros that you will use to type math into
your (roxygen2 or vanilla) documentation. To gain access to these
macros, add the following line to your `DESCRIPTION` file.

``` r
RdMacros: mathjaxr
```

To ensure that users who do not have mathjaxr installed get
documentation with correctly typeset math, add mathjaxr to `Suggests` or
`Imports`. If you add it to `Suggests`, and the user does not have
mathjaxr installed, then MathJax will be sourced from a CDN (via the
Internet) instead of locally from their mathjaxr installation. This
means that math will **not** appear typeset if the user views the HTML
help while offline. To ensure that the user sees typeset math even while
offline, you must add mathjaxr to `Imports`.

Since none of the functions in your R package will actually call
functions from mathjaxr, `R CMD check` and `devtools::check()` will
produce the following note.

``` r
All declared Imports should be used
```

To suppress the note, add the following line to your `NAMESPACE` file.

``` r
import(mathjaxr)
```

## Configuring your documentation files

You must include `\loadmathjax` in all .R (roxygen2) or .Rd (vanilla)
files in which you’d like to type math. If using roxygen2, then
configure your description section like so:

``` r
#' Title text
#' 
#' \loadmathjax
#' Description text
```

(As usual, the `@description` tag has been omitted.) If editing an .Rd
file directly, then this configuration should work:

``` r
\description{
  \loadmathjax
  Description text
}
```

Loading MathJax in the description section (as instructed by the
maintainer) does not prevent you from typing math in other sections,
including `@title`, `@details`, and `@param`. However,
`devtools::document()` will warn

``` r
@tag Use of inline HTML is not currently supported
```

if you have typed inline (not display) math outside of the description
section. I’ve been ignoring these warnings, because the math appears
correctly typeset in the HTML help regardless of the section it is in.

Contradicting the maintainer’s instructions, I have been using this
(roxygen2) configuration without issue:

``` r
#' \loadmathjax
#' Title text
#' 
#' Description text
```

Use at your own risk. I like it because it’s preamble-y.

## Typing math

The workhorse macros are `\mjeqn{}{}` and `\mjdeqn{}{}` for inline math
and display math, respectively. The first argument is the desired
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code *without* dollar ($) delimiters. The second argument
is the desired plain text alternative, which will appear in the plain
text (not HTML) help. If these two arguments happen to be identical
(perhaps you are just typing
![a+b](https://latex.codecogs.com/png.latex?a%2Bb "a+b"), or perhaps you
are fine with raw
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code in the plain text help), then you can use the single
argument macros `\mjseqn{}` and `\mjsdeqn{}` instead.

Here is a simple example from my package oneRous, which implements
several methods for calculating or approximating unity.

``` r
#' Partial sum of negative powers of two
#' 
#' \loadmathjax
#' Calculates
#' \mjdeqn{\sum_{i=1}^{n} \frac{1}{2^i}}{1/2 + ... + 1/(2^n)}
#' for an integer \mjeqn{n \geq 1}{n >= 1}.
#' 
#' Returns unity (1) starting at \mjseqn{n = 54} due to convergence to within
#' machine precision. See [base::.Machine].
#' 
#' @param n An integer scalar. Must be positive.
#' @return A numeric scalar. The value of `sum(cumprod(rep(0.5, n)))`.
#' @export
#' @md
```

A list of supported
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") macros and environments can be found
[here](https://docs.mathjax.org/en/latest/input/tex/macros/index.html).

I have found that macro and environment usage within `\mjeqn{}{}` and
friends can differ from typical usage. For example, I have used the
`align` environment successfully, but I needed to type newlines as `\cr`
instead of `\\`. To get

  
![\\begin{align} x &= r \\cos \\theta \\\\ y &= r \\sin \\theta
\\end{align}](https://latex.codecogs.com/png.latex?%5Cbegin%7Balign%7D%20x%20%26%3D%20r%20%5Ccos%20%5Ctheta%20%5C%5C%20y%20%26%3D%20r%20%5Csin%20%5Ctheta%20%5Cend%7Balign%7D
"\\begin{align} x &= r \\cos \\theta \\\\ y &= r \\sin \\theta \\end{align}")  
in the HTML help, I would write something like
this:

``` r
#' \mjsdeqn{\begin{align} x &= r \cos\theta \cr y &= r \sin\theta \end{align}}
```

The need to type `\cr` instead of `\\` seems to apply more broadly
(i.e., in other multiline environments like `bmatrix`).

## Previewing help without building

Given all of the issues with mathjaxr (see below), it is helpful to be
able to view your package documentation to verify that your math was
typeset, without needing to build and reinstall your package. The
mathjaxr function `preview_rd()` allows you to preview documentation.
For example, after editing `my_function.R`, you can view the updated
HTML help in your browser by running

``` r
devtools::document()
mathjaxr::preview_rd(Rdfile = "my_function.Rd", type = "html")
```

You can view the updated plain text help with `type = "txt"`.

Confusingly, the first argument is just an .Rd file name, not a path to
an .Rd file. You must ensure that your working directory is `man/` or
`man/..`, and specify `Rdfile = "my_function.Rd"` regardless.

## Known issues

Here is an evolving list of issues. Some are mentioned by the maintainer
[on GitHub](https://github.com/wviechtb/mathjaxr), while others I have
discovered. It’s possible that my issues have workarounds that I just
don’t know about.

### Strict inequalities

You must type strict inequalities with spaces (`i < j`, not `i<j`), or
suffer the wrath of the HTML interpreter.

### Literal braces and brackets

Literal braces must be typed as `\lbrace` and `\rbrace` instead of `\{`
and `\}`. Literal brackets must be typed as `\lbrack` and `\rbrack`
instead of `[` and `]`.

### Breaking long lines

Breaking long lines of
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code over multiple lines in my source file caused the
typeset math to vanish. I have addressed this simply by ensuring that
any instance of a mathjaxr macro (e.g., `\mjsdeqn{}`) starts and ends on
the same line.

``` r
#' Good
#' \mjsdeqn{\begin{align} x &= r \cos\theta \cr y &= r \sin\theta \end{align}}
#' 
#' Bad
#' \mjsdeqn{
#'   \begin{align} 
#'     x &= r \cos\theta \cr 
#'     y &= r \sin\theta 
#'   \end{align}
#' }
```

### No newlines in plain text

I have not figured out how to type newlines in plain text, making it
difficult to get something like a system of equations in the plain text
help. Typical escapes like `\n` and `\cr` haven’t worked for me.

``` r
#' Good in HTML help, bad in plain text help
#' \mjdeqn{\begin{align} x - 2y &= 1 \cr 4x + 3y &= 5 \end{align}}{ x - 2y = 1 \cr 4x + 3y = 5}
```

There might be a (cumbersome) workaround involving the `\ifelse` macro,
documented
[here](https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Conditional-text).

### Plain text isn’t protected

When I enable Markdown with `@md`, my plain text gets processed like
Markdown. For example, I wanted `beta*S*I`, but got `beta\emph{S}I` in
the plain text help (literally).

### Macros with subscripts, subscripts with macros (sometimes?)

Some of my equations were typeset fine until I tried to include
something like
![\\widehat{N}\_0](https://latex.codecogs.com/png.latex?%5Cwidehat%7BN%7D_0
"\\widehat{N}_0") (`\widehat{N}_0`),
![\\beta\_\\phi](https://latex.codecogs.com/png.latex?%5Cbeta_%5Cphi
"\\beta_\\phi") (`\beta_\phi`), or
![Z\_\\text{cum}](https://latex.codecogs.com/png.latex?Z_%5Ctext%7Bcum%7D
"Z_\\text{cum}") (`Z_\text{cum}`), at which point raw
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code appeared in the HTML help instead of typeset math.
Curiously, those three terms are typeset fine individually. Including
them as part of a longer line of
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code is what seems to, in some instances but not others,
break the interpreter.

I haven’t deduced exactly what’s going wrong, but I have found that
enclosing all of the
![\\LaTeX{}](https://latex.codecogs.com/png.latex?%5CLaTeX%7B%7D
"\\LaTeX{}") code in `\out{}` solves the issue, except for the fact that
“\\out” is also typeset (literally) before the desired math.

``` r
#' Sometimes good, sometimes bad
#' \mjeqn{Z_\text{cum}}{Z_cum}
#'
#' Good, but "\out" is also typeset
#' \mjeqn{\out{Z_\text{cum}}}{Z_cum}
```

I suspect that there is a workaround without side effects.