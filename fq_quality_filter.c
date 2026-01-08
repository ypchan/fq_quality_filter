// fq_quality_filter.c
// FASTQ/FASTQ.GZ per-read quality fraction filter (multi-threaded workers, ordered output)
// Features requested:
//  1) Log full command line (for reproducibility)
//  2) Time format: "00h 00m 00s"
//  3) Q20/Q30 specialized fast path when rules are exactly {Q>=20, Q>=30} (any fractions)
//  4) stdout output support for piping to pigz, while still supporting built-in single-thread gzip output
//
// Build:
//   gcc -O3 -pipe fq_quality_filter.c -o fq_quality_filter -lz -pthread -Wall -Wextra
// (For portability across different CPUs, avoid -march=native in release builds.)
//
// Examples:
//   # Built-in gzip output (single-thread zlib):
//   ./fq_quality_filter -i in.fastq.gz -o out.fastq.gz --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80
//
//   # Pipe to pigz (recommended for speed):
//   ./fq_quality_filter -i in.fastq.gz -o - --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80 | pigz -p 16 > out.fastq.gz
//
// Notes:
// - Assumes standard 4-line FASTQ records (no wrapped seq/qual lines).
// - Phred auto-detect: samples quality lines and decides 33 vs 64 by min ASCII.
// - Progress: for plain files uses ftello() / file_size; for gz uses gzoffset() (if zlib supports) / gz_size.

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <zlib.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>

#ifndef ZLIB_VERNUM
#define ZLIB_VERNUM 0
#endif

#define LINE_INIT (1<<16)
#define DEFAULT_QUEUE_CAP 4096
#define DEFAULT_SAMPLE_READS 2000
#define MAX_RULES 64

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int ends_with(const char *s, const char *suffix) {
    size_t ls = strlen(s), lf = strlen(suffix);
    if (lf > ls) return 0;
    return strcmp(s + (ls - lf), suffix) == 0;
}

static off_t file_size_bytes(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return (off_t)st.st_size;
}

static void chomp(char *s) {
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == '\n' || s[n-1] == '\r')) {
        s[n-1] = '\0';
        n--;
    }
}

static void fmt_hms(double sec, char *buf, size_t n) {
    if (sec < 0) sec = 0;
    long s = (long)(sec + 0.5);
    long hh = s / 3600;
    long mm = (s % 3600) / 60;
    long ss = s % 60;
    snprintf(buf, n, "%02ldh %02ldm %02lds", hh, mm, ss);
}

static void log_cmdline(int argc, char **argv) {
    fprintf(stderr, "Command: ");
    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        int need_quote = (strchr(a, ' ') || strchr(a, '\t') || strchr(a, '\n'));
        if (need_quote) fputc('"', stderr);
        fputs(a, stderr);
        if (need_quote) fputc('"', stderr);
        if (i + 1 < argc) fputc(' ', stderr);
    }
    fputc('\n', stderr);
}

// ---------- Rules ----------
typedef struct {
    int q;
    double frac;
} QRule;

static int cmp_qrule(const void *a, const void *b) {
    const QRule *x = (const QRule*)a;
    const QRule *y = (const QRule*)b;
    return (x->q < y->q) ? -1 : (x->q > y->q) ? 1 : 0;
}

// ---------- I/O ----------
typedef struct {
    int is_gz;
    gzFile gzf;
    FILE *fp;
    int fd;
    off_t file_size;
    const char *path;
} Reader;

typedef struct {
    int is_gz;
    gzFile gzf;
    FILE *fp;
    const char *path;
} Writer;

static Reader open_reader(const char *path) {
    Reader r;
    memset(&r, 0, sizeof(r));
    r.path = path;

    if (strcmp(path, "-") == 0) {
        r.is_gz = 0;
        r.fp = stdin;
        r.fd = fileno(stdin);
        r.file_size = -1;
        return r;
    }

    r.file_size = file_size_bytes(path);

    if (ends_with(path, ".gz")) {
        r.is_gz = 1;
        r.gzf = gzopen(path, "rb");
        if (!r.gzf) {
            fprintf(stderr, "ERROR: cannot open input gz: %s\n", path);
            exit(1);
        }
        r.fd = open(path, O_RDONLY);
        if (r.fd < 0) r.fd = -1;
    } else {
        r.is_gz = 0;
        r.fp = fopen(path, "rb");
        if (!r.fp) {
            fprintf(stderr, "ERROR: cannot open input: %s\n", path);
            exit(1);
        }
        r.fd = fileno(r.fp);
    }
    return r;
}

static Writer open_writer(const char *path) {
    Writer w;
    memset(&w, 0, sizeof(w));
    w.path = path;

    if (strcmp(path, "-") == 0) {
        w.is_gz = 0;
        w.fp = stdout;
        return w;
    }
    if (ends_with(path, ".gz")) {
        w.is_gz = 1;
        w.gzf = gzopen(path, "wb");
        if (!w.gzf) {
            fprintf(stderr, "ERROR: cannot open output gz: %s\n", path);
            exit(1);
        }
        // Optional: faster compression (bigger files). Uncomment if desired.
        // gzsetparams(w.gzf, 1, Z_DEFAULT_STRATEGY);
    } else {
        w.is_gz = 0;
        w.fp = fopen(path, "wb");
        if (!w.fp) {
            fprintf(stderr, "ERROR: cannot open output: %s\n", path);
            exit(1);
        }
    }
    return w;
}

static void close_reader(Reader *r) {
    if (r->is_gz) gzclose(r->gzf);
    else if (r->fp && r->fp != stdin) fclose(r->fp);
    if (r->fd >= 0 && strcmp(r->path, "-") != 0) close(r->fd);
}

static void close_writer(Writer *w) {
    if (w->is_gz) gzclose(w->gzf);
    else if (w->fp && w->fp != stdout) fclose(w->fp);
}

static void write_str(Writer *w, const char *s) {
    size_t n = strlen(s);
    if (w->is_gz) {
        if ((int)gzwrite(w->gzf, s, (unsigned)n) != (int)n) {
            fprintf(stderr, "ERROR: gzwrite failed\n");
            exit(1);
        }
    } else {
        if (fwrite(s, 1, n, w->fp) != n) {
            fprintf(stderr, "ERROR: write failed\n");
            exit(1);
        }
    }
}

static char* read_line(Reader *r, char **buf, size_t *cap) {
    if (*buf == NULL) {
        *cap = LINE_INIT;
        *buf = (char*)malloc(*cap);
        if (!*buf) { fprintf(stderr, "OOM\n"); exit(1); }
    }

    for (;;) {
        if (r->is_gz) {
            if (gzgets(r->gzf, *buf, (int)(*cap)) == NULL) return NULL;
        } else {
            if (fgets(*buf, (int)(*cap), r->fp) == NULL) return NULL;
        }

        size_t len = strlen(*buf);
        if (len > 0 && (*buf)[len-1] == '\n') return *buf;

        // grow & append remainder
        *cap *= 2;
        *buf = (char*)realloc(*buf, *cap);
        if (!*buf) { fprintf(stderr, "OOM\n"); exit(1); }

        size_t cur = len;
        while (1) {
            char *pos = (*buf) + cur;
            size_t rem = (*cap) - cur;
            if (rem < 2) {
                *cap *= 2;
                *buf = (char*)realloc(*buf, *cap);
                if (!*buf) { fprintf(stderr, "OOM\n"); exit(1); }
                continue;
            }
            if (r->is_gz) {
                if (gzgets(r->gzf, pos, (int)rem) == NULL) break;
            } else {
                if (fgets(pos, (int)rem, r->fp) == NULL) break;
            }
            cur = strlen(*buf);
            if (cur > 0 && (*buf)[cur-1] == '\n') break;
            if (cur + 2 >= *cap) {
                *cap *= 2;
                *buf = (char*)realloc(*buf, *cap);
                if (!*buf) { fprintf(stderr, "OOM\n"); exit(1); }
            }
        }
        return *buf;
    }
}

// ---------- Record / Queue / Ring ----------
typedef struct {
    uint64_t seqno;
    char *h,*s,*p,*q;   // owned
    size_t qlen;
    int keep;
} Record;

static Record* rec_new(uint64_t seqno, const char *h, const char *s, const char *p, const char *q) {
    Record *r = (Record*)calloc(1, sizeof(Record));
    if (!r) { fprintf(stderr, "OOM\n"); exit(1); }
    r->seqno = seqno;
    r->h = strdup(h);
    r->s = strdup(s);
    r->p = strdup(p);
    r->q = strdup(q);
    if (!r->h || !r->s || !r->p || !r->q) { fprintf(stderr, "OOM\n"); exit(1); }
    chomp(r->q);
    r->qlen = strlen(r->q);
    r->keep = 0;
    return r;
}

static void rec_free(Record *r) {
    if (!r) return;
    free(r->h); free(r->s); free(r->p); free(r->q);
    free(r);
}

typedef struct {
    Record **buf;
    int cap, head, tail, count;
    int closed;
    pthread_mutex_t m;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} Queue;

static void q_init(Queue *q, int cap) {
    q->buf = (Record**)calloc((size_t)cap, sizeof(Record*));
    q->cap = cap; q->head = 0; q->tail = 0; q->count = 0; q->closed = 0;
    pthread_mutex_init(&q->m, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
}

static void q_close(Queue *q) {
    pthread_mutex_lock(&q->m);
    q->closed = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_cond_broadcast(&q->not_full);
    pthread_mutex_unlock(&q->m);
}

static void q_destroy(Queue *q) {
    free(q->buf);
    pthread_mutex_destroy(&q->m);
    pthread_cond_destroy(&q->not_empty);
    pthread_cond_destroy(&q->not_full);
}

static void q_push(Queue *q, Record *r) {
    pthread_mutex_lock(&q->m);
    while (q->count == q->cap && !q->closed) pthread_cond_wait(&q->not_full, &q->m);
    if (q->closed) { pthread_mutex_unlock(&q->m); rec_free(r); return; }
    q->buf[q->tail] = r;
    q->tail = (q->tail + 1) % q->cap;
    q->count++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->m);
}

static Record* q_pop(Queue *q) {
    pthread_mutex_lock(&q->m);
    while (q->count == 0 && !q->closed) pthread_cond_wait(&q->not_empty, &q->m);
    if (q->count == 0 && q->closed) { pthread_mutex_unlock(&q->m); return NULL; }
    Record *r = q->buf[q->head];
    q->head = (q->head + 1) % q->cap;
    q->count--;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->m);
    return r;
}

typedef struct {
    uint64_t seqno;
    Record *rec;
    int ready;
} Slot;

typedef struct {
    Slot *slots;
    int cap;
    uint64_t total_in;
    int done_reading;
    pthread_mutex_t m;
    pthread_cond_t cv;
} ResultRing;

static void rr_init(ResultRing *rr, int cap) {
    rr->slots = (Slot*)calloc((size_t)cap, sizeof(Slot));
    rr->cap = cap;
    rr->total_in = 0;
    rr->done_reading = 0;
    pthread_mutex_init(&rr->m, NULL);
    pthread_cond_init(&rr->cv, NULL);
}

static void rr_destroy(ResultRing *rr) {
    free(rr->slots);
    pthread_mutex_destroy(&rr->m);
    pthread_cond_destroy(&rr->cv);
}

// ---------- Phred detection ----------
static int detect_phred(const char *in_path, int sample_reads, int *min_ascii_out, int *max_ascii_out) {
    Reader r = open_reader(in_path);
    char *h=NULL,*s=NULL,*p=NULL,*q=NULL;
    size_t hc=0,sc=0,pc=0,qc=0;

    int min_ascii = 999, max_ascii = -1;
    int n = 0;

    while (n < sample_reads) {
        char *lh = read_line(&r, &h, &hc);
        if (!lh) break;
        char *ls = read_line(&r, &s, &sc);
        char *lp = read_line(&r, &p, &pc);
        char *lq = read_line(&r, &q, &qc);
        if (!ls || !lp || !lq) break;

        chomp(q);
        size_t L = strlen(q);
        for (size_t i=0;i<L;i++) {
            unsigned char c = (unsigned char)q[i];
            if ((int)c < min_ascii) min_ascii = (int)c;
            if ((int)c > max_ascii) max_ascii = (int)c;
        }
        n++;
    }

    free(h); free(s); free(p); free(q);
    close_reader(&r);

    if (min_ascii_out) *min_ascii_out = (min_ascii==999? -1 : min_ascii);
    if (max_ascii_out) *max_ascii_out = max_ascii;

    if (min_ascii == 999) {
        fprintf(stderr, "WARN: cannot sample qualities (empty file?) -> assume phred=33\n");
        return 33;
    }
    if (min_ascii < 59) return 33;
    if (min_ascii >= 64) return 64;

    fprintf(stderr, "WARN: ambiguous quality ASCII_min=%d (59..63). Default phred=33.\n", min_ascii);
    return 33;
}

// ---------- Context ----------
typedef struct {
    // config
    const char *in_path;
    const char *out_path;
    int threads;
    int phred;
    QRule *rules;
    int nrules;
    int sample_reads;

    // fast path for Q20/Q30 rules
    int q20q30_fast;     // 1 if rules are exactly q=20 and q=30
    double q20_frac;
    double q30_frac;

    // io
    Reader in;
    Writer out;

    // pipeline
    Queue q;
    ResultRing rr;

    // counters (protected by rr.m)
    uint64_t reads_read;
    uint64_t reads_done;
    uint64_t reads_kept;

    double start_time;
} Ctx;

// Estimate input position for progress
static off_t input_pos(Reader *r) {
    if (strcmp(r->path, "-") == 0) return -1;

    if (!r->is_gz) {
        return (off_t)ftello(r->fp);
    } else {
#if ZLIB_VERNUM >= 0x1290
        z_off_t off = gzoffset(r->gzf);
        if (off < 0) return -1;
        return (off_t)off;
#else
        return -1;
#endif
    }
}

// ---------- Threads ----------
static void* reader_thread(void *arg) {
    Ctx *ctx = (Ctx*)arg;
    char *h=NULL, *s=NULL, *p=NULL, *q=NULL;
    size_t hc=0, sc=0, pc=0, qc=0;

    uint64_t seqno = 0;
    while (1) {
        char *lh = read_line(&ctx->in, &h, &hc);
        if (!lh) break;
        char *ls = read_line(&ctx->in, &s, &sc);
        char *lp = read_line(&ctx->in, &p, &pc);
        char *lq = read_line(&ctx->in, &q, &qc);
        if (!ls || !lp || !lq) break;

        Record *r = rec_new(seqno, lh, ls, lp, lq);
        q_push(&ctx->q, r);

        pthread_mutex_lock(&ctx->rr.m);
        ctx->reads_read = seqno + 1;
        pthread_cond_broadcast(&ctx->rr.cv);
        pthread_mutex_unlock(&ctx->rr.m);

        seqno++;
    }

    pthread_mutex_lock(&ctx->rr.m);
    ctx->rr.total_in = seqno;
    ctx->rr.done_reading = 1;
    pthread_cond_broadcast(&ctx->rr.cv);
    pthread_mutex_unlock(&ctx->rr.m);

    q_close(&ctx->q);

    free(h); free(s); free(p); free(q);
    return NULL;
}

static inline int eval_rules_generic(const char *qstr, size_t L, const QRule *rules, int nrules, int phred) {
    if (L == 0) return 0;
    int counts[MAX_RULES];
    for (int i=0;i<nrules;i++) counts[i]=0;

    const unsigned char *qs = (const unsigned char*)qstr;
    for (size_t i=0;i<L;i++) {
        int v = (int)qs[i] - phred;
        for (int j=0;j<nrules;j++) {
            if (v >= rules[j].q) counts[j]++;
            else break; // q ascending
        }
    }
    for (int j=0;j<nrules;j++) {
        double frac = (double)counts[j] / (double)L;
        if (frac < rules[j].frac) return 0;
    }
    return 1;
}

static inline int eval_rules_q20q30_fast(const char *qstr, size_t L, int phred, double q20_frac, double q30_frac) {
    if (L == 0) return 0;
    size_t c20 = 0, c30 = 0;
    const unsigned char *qs = (const unsigned char*)qstr;
    for (size_t i=0;i<L;i++) {
        int v = (int)qs[i] - phred;
        if (v >= 20) c20++;
        if (v >= 30) c30++;
    }
    double f20 = (double)c20 / (double)L;
    double f30 = (double)c30 / (double)L;
    return (f20 >= q20_frac && f30 >= q30_frac) ? 1 : 0;
}

static void* worker_thread(void *arg) {
    Ctx *ctx = (Ctx*)arg;

    while (1) {
        Record *r = q_pop(&ctx->q);
        if (!r) break;

        int keep = 0;
        if (ctx->q20q30_fast) {
            keep = eval_rules_q20q30_fast(r->q, r->qlen, ctx->phred, ctx->q20_frac, ctx->q30_frac);
        } else {
            keep = eval_rules_generic(r->q, r->qlen, ctx->rules, ctx->nrules, ctx->phred);
        }
        r->keep = keep;

        pthread_mutex_lock(&ctx->rr.m);
        int idx = (int)(r->seqno % (uint64_t)ctx->rr.cap);
        while (ctx->rr.slots[idx].ready) pthread_cond_wait(&ctx->rr.cv, &ctx->rr.m);

        ctx->rr.slots[idx].seqno = r->seqno;
        ctx->rr.slots[idx].rec = r;
        ctx->rr.slots[idx].ready = 1;
        pthread_cond_broadcast(&ctx->rr.cv);
        pthread_mutex_unlock(&ctx->rr.m);
    }
    return NULL;
}

static void* writer_thread(void *arg) {
    Ctx *ctx = (Ctx*)arg;
    uint64_t expect = 0;

    while (1) {
        pthread_mutex_lock(&ctx->rr.m);

        while (1) {
            int idx = (int)(expect % (uint64_t)ctx->rr.cap);
            if (ctx->rr.slots[idx].ready && ctx->rr.slots[idx].seqno == expect) break;

            if (ctx->rr.done_reading && expect >= ctx->rr.total_in) {
                pthread_mutex_unlock(&ctx->rr.m);
                return NULL;
            }
            pthread_cond_wait(&ctx->rr.cv, &ctx->rr.m);
        }

        int idx = (int)(expect % (uint64_t)ctx->rr.cap);
        Record *r = ctx->rr.slots[idx].rec;

        ctx->rr.slots[idx].ready = 0;
        ctx->rr.slots[idx].rec = NULL;

        ctx->reads_done = expect + 1;
        if (r->keep) ctx->reads_kept++;

        pthread_cond_broadcast(&ctx->rr.cv);
        pthread_mutex_unlock(&ctx->rr.m);

        if (r->keep) {
            write_str(&ctx->out, r->h);
            write_str(&ctx->out, r->s);
            write_str(&ctx->out, r->p);
            write_str(&ctx->out, r->q);
            write_str(&ctx->out, "\n");
        }
        rec_free(r);
        expect++;
    }
}

static void* progress_thread(void *arg) {
    Ctx *ctx = (Ctx*)arg;
    double last = now_sec();

    while (1) {
        usleep(200000); // 0.2s
        double t = now_sec();
        if (t - last < 1.0) continue;
        last = t;

        pthread_mutex_lock(&ctx->rr.m);
        uint64_t done = ctx->reads_done;
        uint64_t kept = ctx->reads_kept;
        int done_reading = ctx->rr.done_reading;
        uint64_t total_in = ctx->rr.total_in;
        pthread_mutex_unlock(&ctx->rr.m);

        double elapsed = t - ctx->start_time;
        double rps = (elapsed > 0) ? (double)done / elapsed : 0.0;

        off_t pos = input_pos(&ctx->in);
        double pct = -1.0;
        if (ctx->in.file_size > 0 && pos >= 0) {
            pct = (double)pos * 100.0 / (double)ctx->in.file_size;
            if (pct > 100.0) pct = 100.0;
        }

        char el_buf[32], eta_buf[32];
        fmt_hms(elapsed, el_buf, sizeof(el_buf));
        strcpy(eta_buf, "NA");

        if (pct >= 0 && elapsed > 1.0 && pct > 0.5 && pct < 99.9) {
            double eta = elapsed * (100.0 - pct) / pct;
            fmt_hms(eta, eta_buf, sizeof(eta_buf));
        }

        if (pct >= 0) {
            fprintf(stderr,
                "\r[progress] %6.2f%%  reads_done=%llu  kept=%llu  r/s=%.1f  elapsed=%s  ETA=%s    ",
                pct,
                (unsigned long long)done,
                (unsigned long long)kept,
                rps, el_buf, eta_buf
            );
        } else {
            fprintf(stderr,
                "\r[progress] reads_done=%llu  kept=%llu  r/s=%.1f  elapsed=%s  ETA=%s    ",
                (unsigned long long)done,
                (unsigned long long)kept,
                rps, el_buf, eta_buf
            );
        }
        fflush(stderr);

        if (done_reading && total_in > 0 && done >= total_in) break;
    }
    fprintf(stderr, "\n");
    return NULL;
}

// ---------- CLI ----------
static void usage(const char *prog) {
    fprintf(stderr,
        "fq_quality_filter\n\n"
        "Usage:\n"
        "  %s -i in.fastq[.gz] -o out.fastq[.gz|-] [--threads 1] \\\n"
        "     -q INT --fraction FLOAT   (repeatable)\n\n"
        "Examples:\n"
        "  %s -i in.fq.gz -o out.fq.gz --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80\n"
        "  %s -i in.fq.gz -o - --threads 8 -q 20 --fraction 0.90 -q 30 --fraction 0.80 | pigz -p 16 > out.fq.gz\n\n"
        "Notes:\n"
        "  - Auto-detect phred 33/64 by sampling quality lines.\n"
        "  - Assumes 4-line FASTQ records.\n",
        prog, prog, prog
    );
}

int main(int argc, char **argv) {
    const char *in_path = NULL;
    const char *out_path = NULL;
    int threads = 1;
    int sample_reads = DEFAULT_SAMPLE_READS;

    QRule rules[MAX_RULES];
    int nrules = 0;

    int have_pending_q = 0;
    int pending_q = -1;

    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--in")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: -i/--in needs value\n"); usage(argv[0]); return 1; }
            in_path = argv[++i];
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: -o/--out needs value\n"); usage(argv[0]); return 1; }
            out_path = argv[++i];
        } else if (!strcmp(argv[i], "--threads")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: --threads needs value\n"); usage(argv[0]); return 1; }
            threads = atoi(argv[++i]);
            if (threads < 1) threads = 1;
        } else if (!strcmp(argv[i], "--sample-reads")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: --sample-reads needs value\n"); usage(argv[0]); return 1; }
            sample_reads = atoi(argv[++i]);
            if (sample_reads < 200) sample_reads = 200;
        } else if (!strcmp(argv[i], "-q")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: -q needs INT\n"); usage(argv[0]); return 1; }
            if (have_pending_q) {
                fprintf(stderr, "ERROR: got -q %d but previous -q %d has no matching --fraction\n", atoi(argv[i+1]), pending_q);
                return 1;
            }
            pending_q = atoi(argv[++i]);
            if (pending_q < 0 || pending_q > 93) {
                fprintf(stderr, "ERROR: -q %d out of range (0..93)\n", pending_q);
                return 1;
            }
            have_pending_q = 1;
        } else if (!strcmp(argv[i], "--fraction")) {
            if (i+1>=argc) { fprintf(stderr,"ERROR: --fraction needs FLOAT\n"); usage(argv[0]); return 1; }
            if (!have_pending_q) {
                fprintf(stderr, "ERROR: --fraction provided without a preceding -q\n");
                return 1;
            }
            double frac = atof(argv[++i]);
            if (frac < 0.0 || frac > 1.0) {
                fprintf(stderr, "ERROR: --fraction %.4f out of range [0,1]\n", frac);
                return 1;
            }
            if (nrules >= MAX_RULES) {
                fprintf(stderr, "ERROR: too many -q/--fraction pairs (max %d)\n", MAX_RULES);
                return 1;
            }
            rules[nrules].q = pending_q;
            rules[nrules].frac = frac;
            nrules++;

            have_pending_q = 0;
            pending_q = -1;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "ERROR: unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (!in_path || !out_path) {
        fprintf(stderr, "ERROR: -i and -o are required.\n");
        usage(argv[0]);
        return 1;
    }
    if (have_pending_q) {
        fprintf(stderr, "ERROR: -q %d has no matching --fraction\n", pending_q);
        return 1;
    }
    if (nrules == 0) {
        fprintf(stderr, "ERROR: must provide at least one pair: -q INT --fraction FLOAT\n");
        return 1;
    }

    // Sort rules; reject duplicates
    qsort(rules, (size_t)nrules, sizeof(QRule), cmp_qrule);
    for (int i=1; i<nrules; i++) {
        if (rules[i].q == rules[i-1].q) {
            fprintf(stderr, "ERROR: duplicate -q %d provided.\n", rules[i].q);
            return 1;
        }
    }

    // Auto-detect phred
    int ascii_min=-1, ascii_max=-1;
    int phred = detect_phred(in_path, sample_reads, &ascii_min, &ascii_max);

    // Decide Q20/Q30 fast path
    int q20q30_fast = 0;
    double q20_frac = 0.0, q30_frac = 0.0;
    if (nrules == 2 && rules[0].q == 20 && rules[1].q == 30) {
        q20q30_fast = 1;
        q20_frac = rules[0].frac;
        q30_frac = rules[1].frac;
    }

    // Log header
    fprintf(stderr, "=== fq_quality_filter start ===\n");
    log_cmdline(argc, argv);
    fprintf(stderr, "Input : %s\n", in_path);
    fprintf(stderr, "Output: %s\n", out_path);
    fprintf(stderr, "Threads(worker): %d\n", threads);
    fprintf(stderr, "Sample reads for phred detection: %d\n", sample_reads);
    fprintf(stderr, "Detected quality ASCII_min=%d ASCII_max=%d -> phred=%d\n", ascii_min, ascii_max, phred);
    fprintf(stderr, "Rules (sorted):\n");
    for (int i=0;i<nrules;i++) fprintf(stderr, "  Q>=%d  fraction>=%.4f\n", rules[i].q, rules[i].frac);
    if (q20q30_fast) fprintf(stderr, "Fast path: enabled (Q20/Q30 specialized)\n");
    else fprintf(stderr, "Fast path: disabled (generic rules)\n");
    fprintf(stderr, "==============================\n");

    // Setup ctx
    Ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.in_path = in_path;
    ctx.out_path = out_path;
    ctx.threads = threads;
    ctx.phred = phred;
    ctx.sample_reads = sample_reads;
    ctx.nrules = nrules;
    ctx.rules = (QRule*)calloc((size_t)nrules, sizeof(QRule));
    if (!ctx.rules) { fprintf(stderr, "OOM\n"); return 1; }
    memcpy(ctx.rules, rules, (size_t)nrules * sizeof(QRule));

    ctx.q20q30_fast = q20q30_fast;
    ctx.q20_frac = q20_frac;
    ctx.q30_frac = q30_frac;

    ctx.in = open_reader(in_path);
    ctx.out = open_writer(out_path);

    int cap = DEFAULT_QUEUE_CAP;
    q_init(&ctx.q, cap);
    rr_init(&ctx.rr, cap);

    ctx.start_time = now_sec();

    // Threads
    pthread_t tr, tw, tp;
    pthread_t *workers = (pthread_t*)calloc((size_t)threads, sizeof(pthread_t));
    if (!workers) { fprintf(stderr, "OOM\n"); return 1; }

    pthread_create(&tr, NULL, reader_thread, &ctx);
    pthread_create(&tw, NULL, writer_thread, &ctx);
    for (int i=0;i<threads;i++) pthread_create(&workers[i], NULL, worker_thread, &ctx);
    pthread_create(&tp, NULL, progress_thread, &ctx);

    pthread_join(tr, NULL);
    for (int i=0;i<threads;i++) pthread_join(workers[i], NULL);

    pthread_mutex_lock(&ctx.rr.m);
    pthread_cond_broadcast(&ctx.rr.cv);
    pthread_mutex_unlock(&ctx.rr.m);

    pthread_join(tw, NULL);
    pthread_join(tp, NULL);

    double elapsed = now_sec() - ctx.start_time;

    char el_buf[32];
    fmt_hms(elapsed, el_buf, sizeof(el_buf));

    fprintf(stderr, "=== fq_quality_filter done ===\n");
    fprintf(stderr, "Total reads   : %llu\n", (unsigned long long)ctx.rr.total_in);
    fprintf(stderr, "Kept reads    : %llu\n", (unsigned long long)ctx.reads_kept);
    fprintf(stderr, "Kept fraction : %.4f\n", ctx.rr.total_in ? (double)ctx.reads_kept/(double)ctx.rr.total_in : 0.0);
    fprintf(stderr, "Elapsed       : %s\n", el_buf);
    fprintf(stderr, "Reads/sec     : %.1f\n", elapsed>0 ? (double)ctx.rr.total_in/elapsed : 0.0);
    fprintf(stderr, "==============================\n");

    // Cleanup
    free(workers);
    free(ctx.rules);
    q_destroy(&ctx.q);
    rr_destroy(&ctx.rr);
    close_reader(&ctx.in);
    close_writer(&ctx.out);

    return 0;
}

