\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{puthesis}[2011/04/29 Princeton University Thesis class Rev. 2.0]

\RequirePackage{setspace}
%color package must be loaded here to avoid conflicts
\RequirePackage[usenames,dvipsnames]{color}

% input definitions
\def\submitted#1{\gdef\@submitted{#1}}
\def\adviser#1{\gdef\@adviser{#1}}
\def\dept#1{\gdef\@dept{#1}}
\def\@deptpref{Department of}
\def\departmentprefix#1{\gdef\@deptpref{#1}}
\def\department#1{\gdef\@dept{#1}}
\long\def\acknowledgements#1{\gdef\@acknowledgements{#1}}
\def\dedication#1{\gdef\@dedication{#1}}
\def\frontnotes#1{\gdef\@frontnotes{#1}}
\long\def\@abstract{\@latex@error{No \noexpand\abstract given}\@ehc}

%generate frontmatter
\newcommand*{\frontmatter}{
  \pagenumbering{roman}
}
\newcommand*{\mainmatter}{\pagenumbering{arabic}}
\newcommand*{\makelot}{}
\newcommand*{\makelof}{}
\newcommand*{\makelos}{}
\newcommand*{\begincmd}{
\doublespacing
\frontmatter\maketitlepage\makecopyrightpage\makeabstract
\makeacknowledgements\makededication\tableofcontents
\clearpage\makelot\clearpage\makelof\clearpage\makelos
\clearpage\makefrontnotes\clearpage\mainmatter }

%the title page
\newcommand{\maketitlepage}{{
  \thispagestyle{empty}
  \sc
  \vspace*{0in}
  \begin{center}
    \LARGE \@title
  \end{center}
  \vspace{.6in}
  \begin{center}
    \@author
  \end{center}
  \vspace{.6in}
  \begin{center}
    A Dissertation \\
    Presented to the Faculty \\
    of Princeton University \\
    in Candidacy for the Degree \\
    of Doctor of Philosophy
  \end{center}
  \vspace{.3in}
  \begin{center}
    Recommended for Acceptance \\
    by the \@deptpref \\
    \@dept \\
    Adviser: \@adviser
  \end{center}
  \vspace{.3in}
  \begin{center}
    \@submitted
  \end{center}
  \clearpage
  }}

%copyright page
\newcommand*{\makecopyrightpage}{
  \thispagestyle{empty}
  \vspace*{0in}
  \begin{center}
    \copyright\ Copyright by \@author, \number\year. \\
    All Rights Reserved
  \end{center}
  \clearpage}
  
%abstract
\newcommand*{\makeabstract}{
  \newpage
  \phantomsection
  \addcontentsline{toc}{section}{Abstract}
  \begin{center}
  \Large \textbf{Abstract}
  \end{center}
  \@abstract
  \clearpage
}

%acknowledgements
\def\makeacknowledgements{
  \ifx\@acknowledgements\undefined
  \else
    \phantomsection
    \addcontentsline{toc}{section}{Acknowledgements}
    \begin{center}
      \Large \textbf{Acknowledgements}
    \end{center}
    \@acknowledgements
    \clearpage
  \fi
  }
  
%dedication
\def\makededication{
  \ifx\@dedication\undefined
  \else
    \vspace*{1in}
    \begin{flushright}
      \@dedication
    \end{flushright}
    \clearpage
  \fi
 }
 
%front notes
\def\makefrontnotes{
  \ifx\@frontnotes\undefined
  \else
    \phantomsection
    \addcontentsline{toc}{section}{Notation}
    \begin{center}
      \Large \textbf{Notation}
    \end{center}
    \@frontnotes
    \clearpage
  \fi
 }

%list of tables
\DeclareOption{lot}{\renewcommand*{\makelot}{
  \phantomsection
  \addcontentsline{toc}{section}{List of Tables}
  \listoftables}}

%list of figures
\DeclareOption{lof}{\renewcommand*{\makelof}{
   \phantomsection
  \addcontentsline{toc}{section}{List of Figures}
  \listoffigures}}

%list of symbols
\DeclareOption{los}{
  \renewcommand*{\makelos}{
    \newcommand*\l@symbol{\@dottedtocline{1}{1.5em}{3.5em}}
    \addtocontents{los}{\protect\addvspace{10\p@}}%
    \chapter*{List of Symbols\@mkboth {LIST OF SYMBOLS}{LIST OF SYMBOLS}}
    \@starttoc{los}
    \addcontentsline{toc}{section}{List of Symbols}
    \newcounter{symbolCounter}
  }
}
%addsymbol command: \addsymbol{symbol}{caption}
\def\addsymbol#1#2{%
         \hypertarget{symbol \arabic{symbolCounter}}{}
	\addtocontents{los}{\protect \contentsline {symbol}{\protect \numberline {#1}{\ignorespaces #2}}{\thepage}{symbol \arabic{symbolCounter}}}
	\stepcounter{symbolCounter}
}

%inherit report class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
\LoadClass{report}

%must load bookmark and geometry packages after the report class
\RequirePackage{bookmark} 
\RequirePackage[left=1.5in,top=1in,right=1in,bottom=1in,footskip=0.333in]{geometry}

%redefine abstract
\long\def\abstract#1{\gdef\@abstract{#1}}

%go!
\AtBeginDocument{\begincmd}
\endinput
%% 
%% End of file `puthesis.cls'.