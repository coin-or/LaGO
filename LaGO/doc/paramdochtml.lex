%option noyywrap

%s comment
%s classknown
%s param

%{
#include <string.h>

int level=0;
char* name=NULL;
char* options=NULL;
char* defaultopt=NULL;
char* descr=NULL;
int printlevel=0;
int classnameprinted=0;
char* classname=NULL;

void printclassname() {
	printf("<H2>%s</H2>", classname);
	free(classname);
	classname=NULL;
}

void printparam() {
	if (level>=printlevel) {
		char* d=descr;
		if (classname) printclassname();
		printf("<TABLE width=100%><TR>\n<TD width=40%><B>%s</B><TD width=40%>",name);
		if (options) printf("%s", options);
		printf("<TD width=20%>");
		if (defaultopt) printf("%s", defaultopt);
		printf("</TR></TABLE>\n");
		while (*d) {
			if (*d=='\n') printf("<BR>");
			printf("%c",*(d++));
		}
	}
	if (name) { free(name); name=NULL; }
	if (options) { free(options); options=NULL; }
	if (defaultopt) { free(defaultopt); defaultopt=NULL; }
}

%}

%%

"/**"	BEGIN(comment);

<comment>@class.*	{ if (classname) free(classname); classname=strdup(yytext+7);  } BEGIN(classknown);

<classknown>@param[ ].*	BEGIN(param); { name=strdup(index(yytext,'@')+7); level=0; }
<param>\%level[ ].*\n { level=atoi(index(yytext, '%')+7); }
<param>\%options[ ].*\n { yytext[yyleng-1]=0; options=strdup(index(yytext, '%')+9); }
<param>\%default[ ].*\n { yytext[yyleng-1]=0; defaultopt=strdup(index(yytext, '%')+9); }
<param>@ { descr=strdup(yytext); descr[yyleng-1]=0; unput('@'); printparam(); } BEGIN(classknown);
<param>"*/"	{ yytext[yyleng-2]=0; descr=strdup(yytext); printparam(); } BEGIN(0);
<param>.|\n { yymore(); }

<*>"*/"	BEGIN(0);

<*>.|\n

%%

main(int argc, char** argv) {
	if (argc>1) printlevel=atoi(argv[1]);
	yylex();
}
