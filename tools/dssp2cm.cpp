/*
** DSSP to regularized CM converter.
**
** Copyright (c) 2007 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int parse_dssp_header(char **seq, char **ss, int **map)
{
	char buf[140];
	/* left uninitialized, these were read below (n_res*n_res for map's
	   calloc) even when the "TOTAL NUMBER OF RESIDUES" line is never found
	   (e.g. empty/non-DSSP input) */
	int n_res = 0, n_chains = 1;
	int found_header = 0;

	while (fgets(buf, sizeof(buf), stdin) != NULL) {
		if (strncmp(buf + 18, "TOTAL NUMBER OF RESIDUES", 24) == 0) {
			sscanf(buf, "%d%d", &n_res, &n_chains);
			n_res += n_chains - 1;
			found_header = 1;
			break;
		}
	}
	if (!found_header)
		fprintf(stderr, "dssp2cm: no \"TOTAL NUMBER OF RESIDUES\" header found; "
				"treating input as an empty structure\n");

	n_res++;			/* extra space make it easier down the road */
	*map = (int *)calloc(n_res * n_res, sizeof(int));
	/* n_res+1: parse_dssp_body indexes seq/ss by the DSSP residue number,
	   which reaches n_res-1, and then writes a terminator one past that. */
	*ss = (char *)calloc(n_res + 1, sizeof(char));
	*seq = (char *)calloc(n_res + 1, sizeof(char));

	while (fgets(buf, sizeof(buf), stdin) != NULL) {
		if (strncmp(buf + 5, "RESIDUE AA STRUCTURE BP1", 24) == 0)
			break;
	}

	return n_res;
}

void parse_dssp_body(int n_res, char *seq, char *ss, int *map)
{
	char buf[140], aa, t, tt = '@', p1, p2, pp1 = ' ', pp2 = ' ';
	/* k = 0: on empty input the loop below never runs and never assigns k,
	   leaving the seq[++k] terminator writing at an indeterminate index. */
	int i, j, k = 0, last_k = 0, cont = -1, c;
	char ibuf[5] = "    ";
	char jbuf[5] = "    ";

	seq[0] = ' ';
	ss[0] = ' ';
	while (fgets(buf, sizeof(buf), stdin) != NULL) {
		if (sscanf(buf, "%5d%*8c%c%*2c%c%*6c%c%c%4c%4c",
			   &k, &aa, &t, &p1, &p2, ibuf, jbuf) != 7)
			continue; //skip chain break
		// k = residue index including chain breaks
		// aa = 1-letter residue ID
		// t = secondary structure type (HEC,BMGIST)
		// p1 = beta strand ladder ID
		// p2 = beta strand ladder ID
		// i = beta-sheet contact partner #1
		// j = beta-sheet contact partner #2
	        if (sscanf(ibuf, "%d", &i) != 1 ||
	            sscanf(jbuf, "%d", &j) != 1)
			continue;
		/* k/i/j come straight from the file and index seq/ss/map, sized off the
		   header's residue count -- a malformed or truncated header lets a body
		   record's index run past the allocation. i/j==0 is the legitimate "no
		   contact" sentinel (see the if(i)/if(j) guards below), so only reject
		   values that would actually be out of range. */
		if (k < 0 || k >= n_res || i < 0 || i >= n_res || j < 0 || j >= n_res) {
			fprintf(stderr, "dssp2cm: residue index %d/%d/%d out of range (0..%d), skipping record\n",
				k, i, j, n_res - 1);
			continue;
		}
		//fprintf(stderr,"%5d        %c  %c      %c%c%4d%4d\n",
		//	   k, aa, t, p1, p2, i, j);

		/* some renaming */
		if (t == 'E') {
			if (i == 0 && j == 0)
				t = ' '; /* strand without interaction partner */
			if (i != 0 && j != 0)
				t = 'M';	/* middle strand */
		}

		if (t == ' ')
			t = 'C';	/* random coils */

		if (aa == '!')
			t = ' '; /* reset type at chain breaks */

		//fprintf(stderr," res %d (%c) ss %c\n",k,aa,t);
		last_k = k; /* the terminator write below must not use a k that got
		               rejected and skipped by the range check above */
		ss[k] = t;
		seq[k] = aa;
		//fprintf(stderr,"%c\n",aa);

		/* detect continuous stretch and determine its length */
		if (t == 'E' || t == 'B' || t == 'M')
			t = 'e';
		else if (t == 'G' || t == 'H' || t == 'I')
			t = 'h';
		else if (t == 'C' || t == 'S')
			t = 'c';
		else if (t == 'T')
			t = 't';

		if (t == tt && (p1 == pp1 || p2 == pp2))
			cont++;
		else
			cont = 0;

		/* sheet contacts */
		if (i) {
			map[k * n_res + i] = map[i * n_res + k] = 1;
			//fprintf(stderr,"beta interaction between %d and %d\n",k,i);
		}

		if (j) {
			map[k * n_res + j] = map[j * n_res + k] = 1;
			//fprintf(stderr,"beta interaction between %d and %d\n",k,j);
		}

		/* helical contacts */
		if (k > 4 && cont >= 4 && t == 'h')
			map[k * n_res + (k - 4)] = map[(k - 4) * n_res + k] = 1;

		if (k > 3 && cont >= 3 && t == 'h')
			map[k * n_res + (k - 3)] = map[(k - 3) * n_res + k] = 1;

		/* second-door neighbors */
		if (k > 2 && cont >= 2) {
			if (t == 'h')
				c = -1;
			else if (t == 'e')
				c = 1;
			else
				c = 0;
			map[k * n_res + (k - 2)] = map[(k - 2) * n_res + k] = c;
		}

		/* next-door neighbors */
		if (k > 1 && cont >= 1) {
			if (t == 'h')
				c = 1;
			else if (t == 'e')
				c = -1;
			else if (t == 'c')
				c = -2;
			else if (t == 't')
				c = 2;
			else
				c = 0;
			map[k * n_res + (k - 1)] = map[(k - 1) * n_res + k] = c;
		}

		/* keep for next round */
		if (aa == '!')
			tt = '@'; /* reset type at chain breaks */
		else
			tt = t;
		pp1 = p1;
		pp2 = p2;

		/* finally its secodary structure */
		map[k * n_res + k] = t;
	}
	seq[++last_k] = '\0';
	ss[last_k] = '\0';
}

void print_map(int n, char *seq, int *map)
{
	int i, j, ch;

	for (i = 1; i < n; putchar('\n'), i++)
		for (j = 1; j < n; j++) {
			if (i == j)
				ch = seq[i];
			else if (map[i * n + j] > 0)
				ch = '+';
			else if (map[i * n + j] == -1)
				ch = '-';
			else
				ch = ' ';
			putchar(ch);
		}
}

void print_special(int n, char *seq, int *map, int skip_chain_break)
{
	int i, j, ch;

	for (i = 1; i < n; i++){

		//skip chain breaks
		if (seq[i] == '!' && skip_chain_break) continue;

		for (j = 1; j < n; j++) {

			//skip chain breaks
			if (seq[j] == '!' && skip_chain_break) continue;

			switch (abs(i - j)) {
			case 0:
				switch (map[i * n + j]) {
				case 'e':
					ch = 'Z';
					break;
				case 'h':
					ch = 'U';
					break;
				default:
					ch = 'O';
				}
				break;
			case 1:
				if (map[i * n + j] > 0 && map[i * n + i] == 'h')
					ch = 'U';
				else if (map[i * n + j] < 0
					 && map[i * n + i] == 'e')
					ch = 'Z';
				else
					ch = 'O';
				break;
			case 2:
				ch = 'O';
				break;
			default:
				if (map[i * n + j] > 0)
					ch = 'X';
				else
					ch = 'O';
			}
			putchar(ch);
		}
		putchar('\n');
	}
}

void print_contacts(int n, char *seq, char *ss, int *map)
{
	int i, j, k;
	char str[11];

	for (i = 1; i < n; i++) {
		str[0] = ss[i];
		str[1] = ' ';
		str[2] = ' ';
		str[3] = (i - 2 > 0 && map[i * n + i - 2]) ? seq[i - 2] : ' ';
		str[4] = (i - 1 > 0 && map[i * n + i - 1]) ? seq[i - 1] : ' ';
		str[5] = seq[i];
		str[6] = (i + 1 < n && map[i * n + i + 1]) ? seq[i + 1] : ' ';
		str[7] = (i + 2 < n && map[i * n + i + 2]) ? seq[i + 2] : ' ';
		str[8] = ' ';
		str[9] = ' ';
		str[10] = '\0';

		/* str[] reserves 2 extra slots each side (8-9 forward, 1-2 backward)
		   beyond the fixed +-1/+-2 neighbours at 3/4/6/7. Without a bound, a
		   residue with >=3 contacts on one side overwrote str[0] or the str[10]
		   terminator, and >=4 wrote outside the array entirely. Capping k stops
		   at the array's real capacity instead -- extra contacts beyond that are
		   dropped, not overflowed. */
		for (j = i + 3, k = 8; j < n && k < 10; j++)
			if (map[i * n + j])
				str[k++] = seq[j];

		for (j = i - 3, k = 2; j > 0 && k >= 0; j--)
			if (map[i * n + j])
				str[k--] = seq[j];

		puts(str);
	}
}

int main(int argc, char *argv[])
{
	int n_res;
	char *ss, *seq;
	int *map;
	int skip_chain_break = 1; //do not add empty rows/columns for chain breaks

	n_res = parse_dssp_header(&seq, &ss, &map);

	parse_dssp_body(n_res, seq, ss, map);
	fprintf(stderr,"seq: %s\n",seq);
	fprintf(stderr," ss: %s\n",ss);

	if (argc > 1 && argv[1][0] == 'c')
		print_contacts(n_res,seq,ss,map);
	else if (argc > 1 && argv[1][0] == 's')
		print_special(n_res, seq, map, skip_chain_break);
	else
		print_map(n_res,seq,map);

	free(map);
	free(ss);
	free(seq);

	return 0;
}
