/*  motif_finder.c — поиск gapped‑мотивов с учётом IUPAC.
    Компиляция:  gcc -std=c99 -O3 -o motif_finder motif_finder.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>

/* --- IUPAC‑коды как битовые маски (A=1, C=2, G=4, T(U)=8) --- */
static const unsigned char IUPAC[256] = {
/*  NUL..  31 */ 0,
 ['A']=1, ['a']=1,
 ['C']=2, ['c']=2,
 ['G']=4, ['g']=4,
 ['T']=8, ['t']=8, ['U']=8, ['u']=8,
 ['R']=1|4, ['r']=1|4,
 ['Y']=2|8, ['y']=2|8,
 ['S']=4|2, ['s']=4|2,
 ['W']=1|8, ['w']=1|8,
 ['K']=4|8, ['k']=4|8,
 ['M']=1|2, ['m']=1|2,
 ['B']=2|4|8, ['b']=2|4|8,
 ['D']=1|4|8, ['d']=1|4|8,
 ['H']=1|2|8, ['h']=1|2|8,
 ['V']=1|2|4, ['v']=1|2|4,
 ['N']=1|2|4|8, ['n']=1|2|4|8,
};

/* расстояние Хэмминга с учётом IUPAC */
static int hamming_iupac(const char *p, const char *t, int m)
{
    int mism = 0;
    for (int i = 0; i < m; ++i)
        if (!(IUPAC[(unsigned char)p[i]] & IUPAC[(unsigned char)t[i]]))
            ++mism;
    return mism;
}

/* левая рукоять: все позиции с ≤ err ошибок */
static void scan_left(const char *seq, int n,
                      const char *left, int mL, int err,
                      int **pos, int *cnt)
{
    *cnt = 0;
    int capacity = 128;
    *pos = malloc(capacity * sizeof(int));
    for (int i = 0; i <= n - mL; ++i) {
        if (hamming_iupac(left, seq + i, mL) <= err) {
            if (*cnt >= capacity) {
                capacity *= 2;
                *pos = realloc(*pos, capacity * sizeof(int));
            }
            (*pos)[(*cnt)++] = i;
        }
    }
}

/* печать FASTA‑заголовка без > и первых пробелов */
static void print_seq_id(const char *header, FILE *out)
{
    const char *p = header;
    while (*p && isspace((unsigned char)*p)) ++p;
    fprintf(out, "%.*s", (int)strcspn(p, " \t\r\n"), p);
}

/* --- CLI --- */
struct opts {
    char  **inputs;   int n_inputs;
    char  *output;
    char  *left;      char *right;
    int    gap_min, gap_max;
    int    err_left, err_right;
    int    after;
};

static void usage(const char *prog)
{
    fprintf(stderr,
"Использование: %s -i in1.fa [in2.fa ..] -o out.tsv \\\n"
"                --left NNN --right NNN --gap-min INT --gap-max INT\n"
"                [--err-left 0] [--err-right 0] [-a AFTER]\n", prog);
    exit(EXIT_FAILURE);
}

static struct opts parse_args(int argc, char **argv)
{
    struct opts o = {0};
    static struct option long_opt[] = {
        {"input",     required_argument, 0, 'i'},
        {"output",    required_argument, 0, 'o'},
        {"left",      required_argument, 0,  1 },
        {"right",     required_argument, 0,  2 },
        {"gap-min",   required_argument, 0,  3 },
        {"gap-max",   required_argument, 0,  4 },
        {"err-left",  required_argument, 0,  5 },
        {"err-right", required_argument, 0,  6 },
        {"after",     required_argument, 0, 'a'},
        {0,0,0,0}
    };
    int c, idx;
    while ((c = getopt_long(argc, argv, "i:o:a:", long_opt, &idx)) != -1) {
        switch (c) {
            case 'i':                         /* список fasta‑файлов */
            o.inputs   = &argv[optind-1]; /* первый уже в optarg  */
            o.n_inputs = 1;               /* ← исправление        */
            /* добавляем всё, что идёт подряд и НЕ начинается с '-' */
            while (optind < argc && argv[optind][0] != '-') {
                ++optind;
                ++o.n_inputs;
            }
            break;
        
        case 'o': o.output = optarg; break;
        case  1 : o.left = optarg; break;
        case  2 : o.right = optarg; break;
        case  3 : o.gap_min = atoi(optarg); break;
        case  4 : o.gap_max = atoi(optarg); break;
        case  5 : o.err_left  = atoi(optarg); break;
        case  6 : o.err_right = atoi(optarg); break;
        case 'a': o.after = atoi(optarg); break;
        default: usage(argv[0]);
        }
    }
    if (!o.inputs || !o.output || !o.left || !o.right)
        usage(argv[0]);
    return o;
}

int main(int argc, char **argv)
{
    struct opts op = parse_args(argc, argv);
    FILE *out = fopen(op.output, "w");
    if (!out) { perror(op.output); return 1; }
    fprintf(out, "seq_id\tstart\tend\tmismatch\tmotif\n");

    /* длины рукояток */
    int mL = strlen(op.left);
    int mR = strlen(op.right);

    /* перебираем все входные FASTA‑файлы */
    for (int f = 0; f < op.n_inputs; ++f) {
        const char *path = op.inputs[f];
        FILE *in = strcmp(path, "-") ? fopen(path, "r") : stdin;
        if (!in) { perror(path); continue; }

        char *line = NULL; size_t len = 0;
        char *seq = NULL;  size_t cap = 0, n = 0;
        char *header = NULL;

        while (getline(&line, &len, in) != -1) {
            if (line[0] == '>') {              /* новый заголовок */
                if (header) goto process_seq;  /* обработать предыдущий */
                header = strdup(line);
                continue;
process_seq:
                /* найти мотивы в накопленной последовательности */
                int *pos = NULL, cnt = 0;
                scan_left(seq, n, op.left, mL, op.err_left, &pos, &cnt);
                for (int p = 0; p < cnt; ++p) {
                    int i = pos[p];
                    int best_end = -1, best_mismR = 999;
                    for (int gap = op.gap_min; gap <= op.gap_max; ++gap) {
                        int j   = i + mL + gap;
                        int end = j + mR - 1;
                        if (end >= n) break;
                        int mismR = hamming_iupac(op.right, seq + j, mR);
                        if (mismR <= op.err_right && mismR < best_mismR) {
                            best_mismR = mismR;
                            best_end   = end;
                            if (!mismR) break;
                        }
                    }
                    if (best_end != -1) {
                        int ext_end = best_end + op.after;
                        if (ext_end >= n) ext_end = n - 1;
                        /* печатаем */
                        print_seq_id(header + 1, out);
                        fprintf(out, "\t%d\t%d\t%d\t",
                            i + 1,
                            ext_end + 1 - op.after,           /* ← правильно */
                            hamming_iupac(op.left, seq + i, mL) + best_mismR);
                        fwrite(seq + i, 1, ext_end - i + 1, out);
                        fputc('\n', out);
                    }
                }
                free(pos); free(seq); seq = NULL; n = cap = 0;
                free(header); header = strdup(line);
            } else {                           /* строка с нуклеотидами */
                size_t l = strcspn(line, "\r\n");
                if (n + l + 1 > cap) {
                    cap = (n + l + 1) * 2;
                    seq = realloc(seq, cap);
                }
                memcpy(seq + n, line, l);
                n += l;
            }
        }
        /* обработать последнюю последовательность */
        if (header && n) goto process_seq;

        free(line); free(seq); free(header);
        if (in != stdin) fclose(in);
    }
    fclose(out);
    return 0;
}
