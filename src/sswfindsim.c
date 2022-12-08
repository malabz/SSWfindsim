#include "sswfindsim.h"
#include "utils.h"
#include <omp.h>

void usage(char name[])
{
	fprintf(stderr, "\n\tUsage: %s [options]\n\n", name);
	fprintf(stderr, "Use Striped-Smith-Waterman method find the similar part in sequence file\n");
	fprintf(stderr, "Available options:\n");
	fprintf(stderr, "\t--center   FILE      center file name (Required)\n");
	fprintf(stderr, "\t--seq      FILE      sequence file name (Required)\n");
	fprintf(stderr, "\t--out      FILE      output file name (Required)\n");
	fprintf(stderr, "\t--matrix   FILE      matrix file name\n");
	fprintf(stderr, "\t--thread   N         threads in this program (default 1)\n");
	fprintf(stderr, "\t--dna                this file is DNA sequence file\n");
	fprintf(stderr, "\t--protein            this file is protein seqeuence file\n");
	fprintf(stderr, "\t--help               print help message\n");
	fprintf(stderr, "\t--version            show program version");
	fprintf(stderr, "Example:\n\t%s --center center.fasta --seq frag.fasta --out frag.txt\n", name);
}

char *c_center, *c_seq, *c_output;
int threads;

kseq_t* centerseq, * seqfile;
int seq_number;
int sequence_type; // -1 is to be determined, 0 is protein, 1 is dna
uint8_t gapopen, gapext;
link_list* head, * nowlist, * nxtlist;
int8_t** sequences_int, *table, *mat, *ref_int;
extern int8_t aa_table[], nt_table[], mat_nt[], blosum50[], blosum62[];
int32_t nmatrix, * readLen, ref_len;

void version()
{
	fprintf(stderr, "sswfindsim Version %d.%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_RELEASE, VER_BUILD);
	exit(0);
}

void get_args(int argc, char* argv[])
{
#if (defined(__linux__) || defined(__APPLE__))
	extern char* optarg;
#elif (defined(_WIN32) || defined(_WIN64))
	extern char* optarg;
#endif
	extern int optind;
	int c;
	threads = 1;
	sequence_type = -1;
	while (1)
	{
		//int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] =
		{
			{"center",    required_argument, 0, 'c'},
			{"out",       required_argument, 0, 'o'},
			{"help",      no_argument,       0, 'h'},
			{"version",   no_argument,       0, 'v'},
			{"seq",       required_argument, 0, 's'},
			{"thread",    required_argument, 0, 't'},
			{"dna",       no_argument,       0, 'd'},
			{"protein",   no_argument,       0, 'p'},
			{0,           0,                 0,  0 }
		};

		c = getopt_long(argc, argv, "c:o:hvs:t:dp", long_options, &option_index);
		if (c == -1) break;
		switch (c)
		{
		case 0:
			fprintf(stderr, "Not supported: option %s", long_options[option_index].name);
			if (optarg) fprintf(stderr, " with arg %s", optarg);
			fprintf(stderr, "\n");
			break;
		case 'h':
			usage(argv[0]);
			version();
			break;
		case 'c':
			c_center = argv[optind - 1];
			break;
		case 'o':
			c_output = argv[optind - 1];
			break;
		case 's':
			c_seq = argv[optind - 1];
			break;
		case 't':
			threads = atoi(argv[optind - 1]);
			if (threads > omp_get_max_threads() || threads <= 0)
			{
				threads = omp_get_max_threads();
				fprintf(stderr, "Warning: Config is not OK. Now use %d thread(s)\n", threads);
			}
			break;
		case 'd':
			sequence_type = 1;
			break;
		case 'p':
			sequence_type = 0;
			break;
		case 'v':
			version();
		}
	}
	if (c_center == NULL || c_output == NULL || c_seq == NULL)
	{
		fprintf(stderr, "Error: Please determine center/output/sequence file. use %s --help for more information.\n", argv[0]);
		exit(EXIT_FAILURE);
	}
}

void init_sequence()
{
	int i;
	fprintf(stderr, "Reading sequences from file %s and %s ... ", c_center, c_seq);
	FILE* f_center = fopen(c_center, "r"), *f_seq = fopen(c_seq, "r");
	if (f_center == NULL || f_seq == NULL) { fprintf(stderr, "\nError: file cannot open. Program will exit\n"); exit(1); }
	centerseq = kseq_init(fileno(f_center));
	seqfile = kseq_init(fileno(f_seq));
	if (kseq_read(centerseq) <= 0) { fprintf(stderr, "\nError: center file %s must have length. Program will exit\n", c_center); exit(1); }
	
	if (sequence_type == -1)
	{
		if (countATGC(centerseq->seq.s) < 0.75) sequence_type = 0;
		else sequence_type = 1;
	}

	if (sequence_type == 0)
	{
		table = aa_table;
		mat = blosum62;
		nmatrix = 24;
	}
	else
	{
		table = nt_table;
		mat = mat_nt;
		nmatrix = 5;
	}

	ref_len = strlen(centerseq->seq.s);
	ref_int = (int8_t*)malloc(sizeof(int8_t) * (size_t)(ref_len + 1));
	if (ref_int == NULL) { fprintf(stderr, "\nError: Can not alloc sequence. Program will exit\n"); exit(1); }
	for (i = 0; i < ref_len; ++i) ref_int[i] = table[(size_t)centerseq->seq.s[i]];

	
	seq_number = 0;
	head = new_node(), nowlist = head;
	while (kseq_read(seqfile) >= 0)
	{
		nxtlist = new_node();
		nowlist->nxt = nxtlist;
		nowlist = nowlist->nxt;
		nowlist->s = (char*)malloc(sizeof(char*) * (strlen(seqfile->seq.s) + 1));
		if(nowlist->s == NULL) { fprintf(stderr, "\nError: Can not alloc sequence. Program will exit\n"); exit(1); }
		strcpy(nowlist->s, seqfile->seq.s);
		++seq_number;
	}
	
	sequences_int = (int8_t**)malloc(sizeof(int8_t*) * (seq_number + 1));
	readLen = (int32_t*)calloc(seq_number + 1, sizeof(int32_t));
	if(sequences_int == NULL || readLen == NULL) { fprintf(stderr, "\nError: Can not alloc sequence matrix. Program will exit\n"); exit(1); }
	int this_id = 0;
	nowlist = head -> nxt;
	while (nowlist)
	{
		readLen[this_id] = strlen(nowlist->s);
		sequences_int[this_id] = (int8_t *) malloc(sizeof(int8_t) * (readLen[this_id] + 1));
		if (sequences_int[this_id] == NULL) { fprintf(stderr, "\nError: Can not alloc sequence matrix. Program will exit\n"); exit(1); }
		for (i = 0; i < readLen[this_id]; ++i) sequences_int[this_id][i] = table[(size_t)nowlist->s[i]];
		++this_id;
		nowlist = nowlist->nxt;
	}
	free_node(head);
	fprintf(stderr, "Done. \n");
}

void swalign()
{
	fprintf(stderr, "Making Stripped Smith-Waterman alignment with %d thread(s) ... \n", threads);
	omp_set_num_threads(threads);
	int i;
	gapopen = 3, gapext = 1;
	s_profile* profile = ssw_init(ref_int, ref_len, mat, nmatrix, 2);
	FILE* f_out = fopen(c_output, "w");
	if (f_out == NULL) { fprintf(stderr, "\nError: Can not open output file %s. Program will exit\n", c_output); exit(1); }
#pragma omp parallel for
	for (i = 0; i < seq_number; ++i)
	{
		s_align* alignment;
		alignment = ssw_align(profile, sequences_int[i], readLen[i], gapopen, gapext, 1, 0, 0, 15);
#pragma omp critical
		alignmentout(alignment, f_out, i);
		align_destroy(alignment);
	}
	fclose(f_out);
	init_destroy(profile);
	fprintf(stderr, "Done. \n");
}

void print_arguments()
{
	fprintf(stderr, "sswfindsim Version %d.%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_RELEASE, VER_BUILD);
	if (sequence_type == 0)
	{
		fprintf(stderr, "Protein sequences, use BLOSUM62 matrix\n");
	}
	else
	{
		fprintf(stderr, "DNA sequences, use simple matrix\n");
	}
	fprintf(stderr, "Stripped Smith Waterman algorithm with gap open penalty = %d, gap extension penalty = %d\n", gapopen, gapext);
}

void freeallfile()
{
	kseq_destroy(centerseq);
	kseq_destroy(seqfile);
}

int main(int argc, char *argv[])
{
	get_args(argc, argv);
	init_sequence();
	swalign();
	print_arguments();
	freeallfile();
	return 0;
}

