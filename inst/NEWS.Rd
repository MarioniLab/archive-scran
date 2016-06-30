\name{scrannews}
\title{scran News}
\encoding{UTF-8}

\section{Version 1.1.4}{\itemize{
\item
Added 'min' argument to testVar().

\item
Transformed correlations to a metric distance in quickCluster().

\item 
Removed normalize() in favour of scater's normalize().

\item 
Allowed setting character strings in isSpike<-.

\item
Allowed compute*Factors() functions to directly return the size factors.

\item 
Added selectorPlot() function for interactive plotting.

\item 
Switched to a group-based weighted correlation for one-way layouts in correlatePairs() and correlateNull(),
and to a correlation of residuals for more complex design matrices.

\item
Added phase assignments to the cyclone() output.

\item
Moved code and subsetting into C++ to improve memory efficiency.

\item
Updated documentation and tests.
}}

\section{Version 1.0.0}{\itemize{
\item
New package scran, for low-level analyses of single-cell RNA sequencing data.
}}

