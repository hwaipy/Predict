/*******************************************************/
/**       PREDICT installation program by KD2BD.      **/
/**   This program is compiled and executed through   **/
/**     the "configure" script in this directory.     **/
/**   It checks for the existence of a soundcard,     **/
/**  creates an appropriate "predict.h" file based    **/
/**  on the directory in which PREDICT was installed, **/
/**  compiles PREDICT and associated utilities, and   **/
/**    sets symbolic links between the executables    **/
/**     generated and the installation directory      **/
/**       specified (/usr/local/bin by default).      **/ 
/**                                                   **/
/**  Created: Oct 1999 -==- Last update: 19-Nov-2010  **/
/*******************************************************/

#include <curses.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/soundcard.h>

char version[10];

void logo()
{
	attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);
	mvprintw(3,18,"                                           ");
	mvprintw(4,18,"         --== PREDICT  v%s ==--         ",version);
	mvprintw(5,18,"   Released by John A. Magliacane, KD2BD   ");
	mvprintw(6,18,"                March 2018                 ");
	mvprintw(7,18,"                                           ");
}

int main (argc,argv)
char argc, *argv[];
{
	int x, cc, dsp=-1;
	char pwd[80], src_path[255], dest_path[255], ans, destination[80];
	FILE *infile, *outfile;

	initscr();
	start_color();
	cbreak();
	scrollok(stdscr,TRUE);

	init_pair(1,COLOR_WHITE,COLOR_BLUE);
	init_pair(2,COLOR_RED,COLOR_BLUE);
	init_pair(3,COLOR_CYAN,COLOR_BLUE);
	init_pair(4,COLOR_GREEN,COLOR_BLUE);
	init_pair(5,COLOR_YELLOW,COLOR_BLUE);
	init_pair(6,COLOR_RED,COLOR_WHITE);

	cc=0;
	getcwd(pwd,79);
	dsp=open("/dev/dsp",O_WRONLY);

	if (dsp!=-1)
		close(dsp);

	infile=fopen(".version","r");
	fscanf(infile,"%s",version);
	fclose(infile);

	bkgdset(COLOR_PAIR(1)|A_BOLD);
	clear();
	refresh();
	logo();
	attrset(COLOR_PAIR(5)|A_BOLD);

	mvprintw(10,2,"PREDICT is a satellite tracking and orbital prediction program written for\n");
	printw("  Linux and similar operating systems by John A. Magliacane, KD2BD.\n");
	printw("  PREDICT is free software.  You can redistribute it and/or modify it under\n");
	printw("  the terms of the GNU General Public License as published by the Free\n");
	printw("  Software Foundation, either version 2 of the License or any later version.\n\n");
	printw("  PREDICT is distributed in the hope that it will useful, but WITHOUT ANY\n");
	printw("  WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS\n");
	printw("  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more\n");
	printw("  details.\n\n");

	attrset(COLOR_PAIR(1)|A_BOLD);
	printw("  Do you accept these conditions and wish to install this software? [Y/N] ");
	refresh();

	do
	{
		ans=getch();

	} while (ans!='n' && ans!='N' && ans!='y' && ans!='Y' && ans!=27);


	if (ans=='y' || ans=='Y')
	{
		clear();
		curs_set(0);
		logo();
		attrset(COLOR_PAIR(4)|A_BOLD);
		printw("\n\n\n\n  PREDICT appears to be installed under %s/\n",pwd);

		if (dsp!=-1)
			printw("  An audio device was found at /dev/dsp");
		else
			printw("  No soundcard was found in your system... Bummer!");

		/* Write predict.h */

		outfile=fopen("predict.h","w");

		fprintf(outfile,"/* This file was generated by the installer program */\n\n");

		fprintf(outfile,"char *predictpath={\"%s/\"}, ",pwd);

		if (dsp==-1)
			fprintf(outfile, "soundcard=0,");
		else
			fprintf(outfile, "soundcard=1,");

		fprintf(outfile, " *version={\"%s\"};\n",version);

		fclose(outfile);

		printw("\n  predict.h was successfully created!\n");
		attrset(COLOR_PAIR(3)|A_BOLD);
		printw("\n  Now compiling PREDICT.  (This may take a while...)\n");
		refresh();

		/* Compile PREDICT... */

		cc=system("cc -Wall -O3 -s -fomit-frame-pointer -lm -lncurses -pthread predict.c -o predict");

		/* Create vocalizer.h */

		outfile=fopen("vocalizer/vocalizer.h","w");

		fprintf(outfile,"/* This file was generated by the installer program */\n\n");
		fprintf(outfile,"char *path={\"%s/vocalizer/\"};\n",pwd);

		fclose(outfile);

		if (dsp!=-1 && cc==0)
		{
			attrset(COLOR_PAIR(3)|A_BOLD);
			printw("  Compiling vocalizer...\n");
			refresh();
			system("cc -Wall -O3 -s -fomit-frame-pointer vocalizer/vocalizer.c -o vocalizer/vocalizer");
		}

		/* Now install the programs and man page by creating symlinks
		   between the newly created executables and the destination
		   directory.  The default destination directory for the
	   	   executables is /usr/local/bin.  This default may be
		   overridden by specifying a different path as an argument
		   to this program (ie: installer /usr/bin).  Normally this
		   is passed along from the "configure" script. */
	
		if (argc==2)
			strncpy(destination,argv[1],78);
		else	
			strncpy(destination,"/usr/local/bin/\0",16);

		/* Ensure a trailing '/' is
		   present in "destination". */

		x=strlen(destination);

		if (destination[x-1]!='/' && x!=0)
		{
			destination[x]='/';
			destination[x+1]=0;
		}

		if (cc==0)
		{
			attrset(COLOR_PAIR(3)|A_BOLD);
			printw("  Linking PREDICT binaries to %s\n\n",destination);
			sprintf(dest_path,"%spredict",destination);
			unlink(dest_path);
			sprintf(src_path,"%s/predict",pwd);
			symlink(src_path,dest_path);
			sprintf(dest_path,"%sxpredict",destination);
			unlink(dest_path);
			sprintf(src_path,"%s/xpredict",pwd);
			symlink(src_path,dest_path);
			unlink("/usr/local/man/man1/predict.1");
			sprintf(dest_path,"%s/docs/man/predict.1",pwd);
			symlink(dest_path,"/usr/local/man/man1/predict.1");

			attrset(COLOR_PAIR(5)|A_BOLD);
			printw("  Don't forget to check out the new graphical\n");
			printw("  client applications under the 'clients' directory!");
			attrset(COLOR_PAIR(1)|A_BOLD);
			printw("\n\n  Done!  Visit http://www.qsl.net/kd2bd/predict.html for the latest news!");
		}

		else
		{
			attrset(COLOR_PAIR(2)|A_BOLD);
			printw("  *** Compilation failed.  Program not installed.  :-(");
			beep();
		}
	}

	refresh();
	unlink("installer");
	curs_set(1);	
	refresh();
	endwin();
	exit(0);
}
