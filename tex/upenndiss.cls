% U. Penn dissertation class.  Version 0.2, 5/24/2000
% Reimplementation based on penndiss.sty (a local adaptation of suthesis.sty)
% Alexis Dimitriadis 12/28/99 (alexis@ling.upenn.edu)
%
% A sample input file is appended below, after the command \endinput
%
% Use: \documentclass{upenndiss}
% Options: 
%   Most report class options; no additional ones.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{upenndiss}[2000/5/24 U. Penn Dissertation Class for LaTeX2e]

% The default point size is 12pt, not 10pt.
% Otherwise the mandated margins result in lines that are much too long
% You can specify another size with, e.g., the [11pt] class option.
\def\@sizeopt{12pt}
\DeclareOption{10pt}{\def\@sizeopt{10pt}}
\DeclareOption{11pt}{\def\@sizeopt{11pt}}
 
% Ignore [(no)titlepage] option (should we warn?)
\DeclareOption{titlepage}{\relax}
\DeclareOption{notitlepage}{\relax}
% I used to pass [notitlepage] to the report class, but its effects seem
% to have been voided by the reimplementation of the front matter.
% (Restore if it causes any problems)
%\PassOptionsToClass{notitlepage}{report}

% Pass all other options to the parent class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions 

% Invoke the parent class with appropriate font size (default 12pt).
% Options given to \PassOptionsToClass above are implicitly passed.
\LoadClass[\@sizeopt]{report}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           START OF CUSTOMIZABLE SECTION
%
% Redefining or overriding these settings should allow adaptation of
% this style file to other dissertation specifications.
%

% Dissertation margins (based on U.S. Letter paper):

% Margins: 1.5in inside margin, 1in outside, top, bottom
\textwidth=6in
\oddsidemargin=0.5in \evensidemargin=0.5in  
\topmargin=0.0in
\textheight=9in  % 9in must include page numbers
\advance\textheight by -1\footskip \relax

\usepackage{setspace}
\doublespacing

%%%% Content of text labels; can be overridden from document preamble  %%%%
\def\acknowledgementsname{Acknowledgements}
\def\prefacename{Preface}

% Position titles given with the signatures on the title page
\def\supervisorname{Supervisor of Dissertation}
\def\approvedname{Graduate Group Chair}
\def\committeemembername{Committee Member}

% Select one in your document: \copyrighttrue or \copyrightfalse
% Controls whether to create a copyright page.
\newif\ifcopyright \copyrighttrue 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         Definitions of front matter contents

% Front matter elements are defined in the preamble and assembled in
% the proper order by \FrontMatter.

% The following commands should be used in the preamble to provide the
% appropriate text or \input commands.  Some are optional; if not sure,
% you can always just try leaving them out. The style will complain if they
% are obligatory!

\def\department#1{\gdef\@department{#1}}
\def\supervisor#1{\gdef\@supervisor{#1}} % Separate joint supervisors with \and
\let\superviser=\supervisor 
\let\advisor=\supervisor
\let\supervisors=\supervisor
\let\supervisers=\supervisor

% The Penn linguistics department lists committee names,
% so here they are; (optional command)
\def\committee#1{\gdef\@committee{#1}} % Separate committee members with \and

\def\approvedby#1{\gdef\@approvedby{#1}}
\let\gradchair=\approvedby

\def\copyrightyear#1{\gdef\@copyrightyear{#1}}  
\copyrightyear{\number\year}         % may be overridden in the document


% Except maybe for \dedication, the following are typically too long to
% appear as in-line commands, so they are not defined as \long
% to encourage the use of separate files.
% Example: \dedication{blah, blah} OR \dedicationfile{Dedication}

% The dedication is vertically centered on the page, and otherwise
% formatted by the author
\long\def\dedication#1{%
    \long\def\@dedication{{\null\vfil\parindent=0pt\relax #1\par}}}
\def\acknowledgements#1{%
    \def\@acknowledgements{\PreSection{\acknowledgementsname} #1}}
\def\abstract#1{\def\@abstract{#1}}
\def\preface#1{\def\@preface{\PreSection{\prefacename} #1}}

% Contrive to \include a file without starting a new page (by temporarily
% disabling \clearpage).  A new page is still begun at the end of the include.
\def\@includehere#1{\bgroup\let\clearpage=\egroup\include{#1}}

\def\dedicationfile#1{\dedication{\@includehere{#1}}}
\def\abstractfile#1{\abstract{\@includehere{#1}}}
\def\acknowledgementsfile#1{\acknowledgements{\@includehere{#1}}}
\def\prefacefile#1{\preface{\@includehere{#1}}}

% The format and order of the preliminary pages.
% Optional pages are only created if requested.
% List of Figures and List of Tables pages are only created if there are
% figures or tables. 
% Dedication, acknowledgements and abstract macros are cleared after use
% to save space.
% Syntax: 
%    \FrontMatter 
% or 
%    \FrontMatter* 
%    <additional preliminary sections, e.g., \include{Abbreviations}> 
%    \MainMatter
%
\def\FrontMatter{%
  \pagenumbering{roman}\pagestyle{plain}%
  \@titlepage 
  \ifcopyright \@copyrightpage \fi 
  \@ifundefined{@dedication} \relax
      {\@dedication \newpage \let\@dedication=\@undefined}%
  \@ifundefined{@acknowledgements} \relax
      {\@acknowledgements \newpage \let\@acknowledgements=\@undefined}%
  % The abstract and TOC are obligatory
  \@abstractpage \newpage \let\@abstract=\@undefined
  \tableofcontents
  \@ifundefined{pd@list@tables} \relax 
      \listoftables 
%c:   Now done automatically
%c    {\addcontentsline{toc}{chapter}{\listtablename}\listoftables}%
%  \@ifundefined{pd@list@figures} \relax 
      \listoffigures
%c    {\addcontentsline{toc}{chapter}{\listfigurename}\listoffigures}%
  \@ifundefined{@preface} \relax 
      {\@preface \newpage \let\@preface=\@undefined}%
  \@ifstar{}{\MainMatter}}


% Start the regularly-numbered part of the dissertation.  Must be given
% explicitly if \FrontMatter* was used.
%
\def\MainMatter{\newpage\pagenumbering{arabic}}


% \@titlepage is modular so that its parts can be redefined as necessary.
% The main macro sets the title and names; \@signatureblock creates the
% signature lines.

\def\@titlepage{%
  \thispagestyle{empty}\enlargethispage{\the\footskip}%
  \null\vskip.1in%
  \begin{center}
%	{\setstretch{2.5} \uppercase\expandafter{\@title}\par }%
	{\setstretch{2.5} \MakeUppercase{\@title}\par }%
	\vskip.3in
	\@author
	\vskip .3in
	A DISSERTATION \\[.1in]
	in \\[.1in]
	\@department
	\vfill
Presented to the Faculties of the University of Pennsylvania in Partial \\
Fulfillment of the Requirements for the Degree of Doctor of Philosophy
\\[0.3in]
	\@copyrightyear
  \end{center}
  \bgroup \parindent=0pt \par \@signatureblock \egroup
  \newpage}
% (There is no point to making ``University of Pennsylvania'' a macro: other
% schools are sure to have a different format for the title page anyway).


% Format the signatures at the bottom of the title page.

\def\@signatureblock{\vskip 0pt plus 0.5 fill
  % Count the total of signatures needed (total names in \@supervisor,
  % \@committee, and \@approvedby)
  \global\@tempcnta=0 
  \@count@names\@tempcnta\@supervisor
  \@count@names\@tempcnta\@committee
  \global\advance\@tempcnta by 1 \relax % for the grad chair or approver
  % Choose mode: one or two columns of signatures (laid out row by row)
  % Then generate the signatures in the desired order 
  \ifnum\@tempcnta < 4
    \let\@eachsig=\@widesig
    \@siglist\@supervisor\supervisorname
    \@siglist\@committee\committeemembername
    \@siglist\@approvedby\approvedname
  \else 
    \let\@eachsig=\@narrowsig
    \@siglist\@supervisor\supervisorname
    \@siglist\@approvedby\approvedname
    \@siglist\@committee\committeemembername
  \fi }

% If the default signature layout is inappropriate, just override it:
% Redefine \@signatureblock and lay out the names by hand, either using
% \@widesig and/or \@narrowsig, or just drawing lines in your favorite way.
% Example: 
% \renewcommand\@signatureblock{%
%   \@widesig{John Doe \\ Supervisor of Dissertation}%
%   \@widesig{Approved by: Jane Blow, Minister of Truth}}


% Signature format intended for signatures laid out one per line
\def\@widesig#1{\parbox[t]{12cm}{\parindent=0pt \vskip 0.65in
  \par \hrule width 7cm \vskip .1in #1}\\}

% This is narrow enough that it should come out two per line unless
% followed by \\
\def\@narrowsig#1{\parbox[t]{7cm}{\parindent=0pt \vskip 0.65in
  \raggedright
  \par \hrule width 7cm \vskip .1in #1}\hfill}


% Do a sig line for each name in the list #1, with title #2
% The output is formatted in the format Name, Title (you can change it here)
% Example: \@siglist{\@committee}{Committee Member}
\def\@siglist#1#2{\def\@aloopuse##1{\@eachsig{##1, \penalty-1 #2}}\@alooprun#1}

% Alternative format: put title below the name, with no comma
% \def\@siglist#1#2{\def\@aloopuse##1{\@eachsig{##1\\ #2}}\@alooprun#1}

% The format of the copyright page
\def\@copyrightpage{%
  \thispagestyle{empty}%
  \null\vfill
  \begin{center}
    COPYRIGHT \\[0.3in]
    \@author \\[0.3in]
    \@copyrightyear
  \end{center}
  \vfill\newpage
}

% The format of the abstract heading
\def\@abstractpage{\addcontentsline{toc}{chapter}{Abstract}%
  \begin{center}
    {\large\bf ABSTRACT} \\[6pt]
    \MakeUppercase{\@title} \\[4pt]
%    Author: \@author \\[4pt]
    \@author \\[4pt]
    \def\and{and }%
    \global\@tempcnta=0 \@count@names\@tempcnta\@supervisor
    \ifnum \@tempcnta > 1 Supervisors: \else Supervisor: \fi \@supervisor
  \end{center}
  \vskip 10pt
  \par\noindent \@abstract }



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  END OF CUSTOMIZABLE SECTION
%
%  LaTeX magic beyond this point-- Use extreme care if you have to
%  make changes!
%

% Looping macros for name lists.
% The elements of \arglist must be separated by <space>\and
% Syntax: \@alooprun\arglist
%
\def\@alooprun#1{\ifx #1\undefined \else 
  \bgroup \let\@aloopend=\egroup 
  \def\and{\@ifnextchar\@aloopend \relax \@loopa}%
  \def\@loopa##1 \and{\@aloopuse{##1}\and}%
  \expandafter\and #1 \and\@aloopend \fi}

% Define \@aloopuse to do whatever we want with each partial argument
% Example: \def\@aloopuse#1{\@eachsig{#1, Committee Member}}

% Count the names in the \and-separated list #2, and ADD the result into 
% the counter #1
% Syntax: \@count@names\counter\list
\def\@count@names#1#2{%
  \def\@aloopuse##1{\global\advance#1 by 1\relax}% count and discard argument
  \@alooprun#2}


% Make a note of whether there were figures and tables
\let\up@svd@caption=\caption
\def\caption
  {\expandafter\gdef\csname pd@have\@captype s\endcsname{1} \up@svd@caption}

% The following code appends; I need to prepend or the implicit argument of
% \caption (the name of the figure or table) is lost!
%% (no \protect needed when defining directly)
%\g@addto@macro\caption
%  {\protect\expandafter\gdef\csname pd@have\@captype s\endcsname{1}}

% Save a *different* flag in .aux (so it can turn off if the tables go away).
\AtEndDocument{
  \@ifundefined{pd@havetables} \relax
     {\immediate\write\@mainaux{\string\gdef\string\pd@list@tables{1}}}
  \@ifundefined{pd@havefigures} \relax
     {\immediate\write\@mainaux{\string\gdef\string\pd@list@figures{1}}}}


% Start a section in the front matter, adding it to the TOC
% Syntax: \PreSection[TOC title]{Actual title}
%
\def\PreSection{\chapter*}
% (This is trivial thanks to the following modification to default LaTeX
% behavior: 


% Modify handler for \chapter* so that the title is automatically added to the
% toc.  Lines obviated by this in earlier code are commented out with %c.
% Let me know if this feature causes problems; it can be removed if
% necessary.

% This also extends the \chapter[TOCtitle]{actual title} syntax to \chapter*.
%
\let\orig@@schapter=\@schapter
% If no optional argument in [ ] follows, \@dblarg copies the obligatory one
\def\@schapter{\@dblarg{\up@schapter}}
\def\up@schapter[#1]#2{\addcontentsline{toc}{chapter}{#1}\orig@@schapter{#2}}

% If you need to restore the default behavior of \chapter*, just say:
%     \let\@schapter=\orig@@schapter
% You should almost certainly not do this until after the front matter!


\endinput







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	   END OF THE UPENNDISS.CLS STYLE

Copy everything below this line into a separate file to create a
sample input file for upenndiss.cls:


%-------------(cut here)-----------------------------------------------



% Sample document for  UPENNDISS.CLS, a document class for University of
% Pennsylvania Ph.D. dissertations (but especially tailored to linguistics).
%
% Alexis Dimitriadis 12/28/99 (alexis@ling.upenn.edu)
%


\documentclass{upenndiss}
% The default font size is 12pt, but can be changed in the normal way:
% \documentclass[11pt]{upenndiss}
%
% The Penn guidelines allow smaller sizes, but be warned that the lines
% are too long for anything smaller than 12pt to be easily readable (see any
% introduction to designing text layouts).

% Add packages and definitions you want to use here:
\usepackage{times}


%%%%%%%%%%%%%%%%%%%  The title pages  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Use \protect\\ if you need a manual line break
\title{A Really Clever Title: \protect\\ The Second Line of Your Title}

% Fill in here as appropriate.  Separate multiple names with \and  
\author{Your Name Here}
\supervisor{Anthony Kroch}
\gradchair{Donald A. Ringe}
\committee{One Professor \and Another Professor}

% Change to \copyrighttrue if you want to generate a copyright page.
\copyrightfalse  
% \copyrightyear{2000}  % Defaults to the current year if undefined.

\department{Linguistics}


%%%%%%%%%%%%%%%%%%%  Now for the front matter %%%%%%%%%%%%%%%%%%%%%%%%%

% Front matter elements are required to appear in a specific order.  In this
% section you define the parts of your front matter, and the style will
% automatically put them in the proper order.


% Declare the sections that you use.  Only Abstract is obligatory.  
%
\abstractfile{Abstract}
\acknowledgementsfile{Acknow}
% \dedicationfile{dedication}
% \prefacefile{preface}

% ALTERNATELY, if you have a really short Dedication (no more than one
% paragraph) you can place it in the main file, with the command
% \dedication{...}, INSTEAD of using \dedicationfile.
\dedication{\center
This example is dedicated to all those users \\
who hate to read manuals.
}

\begin{document}
% The following command now generates all your starting pages:
\FrontMatter

\chapter{Introduction}

This is a sample document for \textsf{upenndiss.cls,} a style class for
University of Pennsylvania Ph.D. dissertations (but especially tailored to the
needs of linguistics dissertations).  It meets the style requirements given in
the University of Pennsylvania Doctoral Dissertation Manual
(\emph{http://www.upenn.edu/VPGE/DissManual.html}) as of December, 2000. 

You can use this document as a template: Delete or ignore the comments,
and adapt the commands to your own needs.  Normally, there would be no
real text in this file: you should place each chapter in a separate file (with
its name ending in \textsf{.tex}, and use a \verb|\include| command to include
it, as below.

Page numbering is automatically changed from roman to regular numerals after
the front matter.  What if you want to add something to the front matter, such
as a list of abbreviations?  Then do the following:

\begin{verbatim}
\FrontMatter*
\include{ListOfSymbols}
\MainMatter
\include{Introduction}
etc.
\end{verbatim}

More ways to customize the format are described in comments in the class file,
\textsf{upenndiss.cls.}  Look there if you need to modify its behavior.

The class file only deals with the formatting requirements given in the
Doctoral Dissertation Manual.  You must make your own arrangements for
bibliographic citations, numbered examples, etc. (You could look at
\emph{http://www.ling.upenn.edu/advice/latex/} for some pointers).


\include{Chapter2}
\include{Chapter3}
\include{Conclusions}


% And the bibliography

\bibliographystyle{linquiry2}
\bibliography{myreferences}

\end{document}


