/* huff.c
 * Huffman compression/decompression using:
 * - Linear data structure: singly linked list (priority queue)
 * - Non-linear data structure: binary tree (Huffman tree)
 *
 * Format:
 *  [4 bytes]  Magic "HUF1"
 *  [8 bytes]  Original file size (uint64_t, little-endian)
 *  [256*4]    Frequency table of bytes (uint32_t, little-endian)
 *  [payload]  Huffman-encoded bitstream (MSB-first within bytes)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

/* -------------------- Utilities -------------------- */

static void die(const char *msg) {
    perror(msg);
    exit(EXIT_FAILURE);
}

static void die_msg(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

static size_t fread_or_die(void *ptr, size_t size, size_t nmemb, FILE *f, const char *ctx) {
    size_t got = fread(ptr, size, nmemb, f);
    if (got != nmemb) {
        if (ferror(f)) die(ctx);
        /* EOF allowed in some contexts; caller should check returned size */
    }
    return got;
}

static void fwrite_or_die(const void *ptr, size_t size, size_t nmemb, FILE *f, const char *ctx) {
    if (fwrite(ptr, size, nmemb, f) != nmemb) die(ctx);
}

/* -------------------- Data Structures -------------------- */

/* Non-linear structure: Huffman Tree Node */
typedef struct Node {
    uint32_t freq;
    int is_leaf;
    uint16_t symbol;       /* 0..255 if leaf */
    struct Node *left;
    struct Node *right;
} Node;

/* Linear structure: singly linked list node for priority queue */
typedef struct ListNode {
    Node *tree;
    struct ListNode *next;
} ListNode;

/* Create a new tree node */
static Node* make_node(uint32_t freq, int is_leaf, uint16_t symbol, Node *l, Node *r) {
    Node *n = (Node*)malloc(sizeof(Node));
    if (!n) die("malloc");
    n->freq = freq;
    n->is_leaf = is_leaf;
    n->symbol = symbol;
    n->left = l;
    n->right = r;
    return n;
}

/* Insert tree node into sorted linked list by ascending freq */
static void list_insert_sorted(ListNode **head, Node *tree) {
    ListNode *ln = (ListNode*)malloc(sizeof(ListNode));
    if (!ln) die("malloc");
    ln->tree = tree;
    ln->next = NULL;

    if (!*head || (*head)->tree->freq > tree->freq) {
        ln->next = *head;
        *head = ln;
        return;
    }
    ListNode *cur = *head;
    while (cur->next && cur->next->tree->freq <= tree->freq) {
        cur = cur->next;
    }
    ln->next = cur->next;
    cur->next = ln;
}

/* Pop the smallest-frequency tree from list */
static Node* list_pop_front(ListNode **head) {
    if (!*head) return NULL;
    ListNode *ln = *head;
    Node *t = ln->tree;
    *head = ln->next;
    free(ln);
    return t;
}

/* Free entire Huffman tree */
static void free_tree(Node *root) {
    if (!root) return;
    free_tree(root->left);
    free_tree(root->right);
    free(root);
}

/* -------------------- Bit I/O -------------------- */

typedef struct {
    FILE *f;
    uint8_t buf;
    int bits_filled; /* number of bits already in buf [0..7] */
} BitWriter;

static void bw_init(BitWriter *bw, FILE *f) {
    bw->f = f;
    bw->buf = 0;
    bw->bits_filled = 0;
}

/* Write one bit (0/1), MSB-first in output bytes */
static void bw_write_bit(BitWriter *bw, int bit) {
    bw->buf = (bw->buf << 1) | (bit ? 1 : 0);
    bw->bits_filled++;
    if (bw->bits_filled == 8) {
        fwrite_or_die(&bw->buf, 1, 1, bw->f, "fwrite(bit)");
        bw->buf = 0;
        bw->bits_filled = 0;
    }
}

static void bw_write_bits(BitWriter *bw, uint32_t code, int len) {
    /* code's most significant (len) bits are written first */
    for (int i = len - 1; i >= 0; --i) {
        int bit = (code >> i) & 1u;
        bw_write_bit(bw, bit);
    }
}

static void bw_flush(BitWriter *bw) {
    if (bw->bits_filled > 0) {
        /* pad remaining with zeros in LSB */
        bw->buf <<= (8 - bw->bits_filled);
        fwrite_or_die(&bw->buf, 1, 1, bw->f, "fwrite(flush)");
        bw->buf = 0;
        bw->bits_filled = 0;
    }
}

typedef struct {
    FILE *f;
    uint8_t buf;
    int bits_left; /* number of unread bits in buf [0..8] */
} BitReader;

static void br_init(BitReader *br, FILE *f) {
    br->f = f;
    br->buf = 0;
    br->bits_left = 0;
}

/* Returns 0/1 for a bit, and sets *ok=0 on EOF without a full byte */
static int br_read_bit(BitReader *br, int *ok) {
    if (br->bits_left == 0) {
        size_t got = fread(&br->buf, 1, 1, br->f);
        if (got != 1) {
            *ok = 0;
            return 0;
        }
        br->bits_left = 8;
    }
    int bit = (br->buf >> 7) & 1;
    br->buf <<= 1;
    br->bits_left--;
    *ok = 1;
    return bit;
}

/* -------------------- Huffman Core -------------------- */

/* Build Huffman tree from frequency table. Uses linked list as min-PQ. */
static Node* build_huffman_tree(const uint32_t freq[256]) {
    ListNode *pq = NULL;
    int symbols = 0;
    for (int i = 0; i < 256; ++i) {
        if (freq[i] > 0) {
            Node *leaf = make_node(freq[i], 1, (uint16_t)i, NULL, NULL);
            list_insert_sorted(&pq, leaf);
            symbols++;
        }
    }
    if (symbols == 0) {
        return NULL; /* empty file: no tree */
    }
    if (symbols == 1) {
        /* Special: single symbol; create a parent to ensure at least one bit */
        Node *only = list_pop_front(&pq);
        Node *root = make_node(only->freq, 0, 0, only, NULL);
        return root;
    }
    while (pq && pq->next) {
        Node *a = list_pop_front(&pq);
        Node *b = list_pop_front(&pq);
        Node *parent = make_node(a->freq + b->freq, 0, 0, a, b);
        list_insert_sorted(&pq, parent);
    }
    Node *root = list_pop_front(&pq);
    return root;
}

/* Code table: for each byte, a (code, length) pair */
typedef struct {
    uint32_t code;
    uint8_t  len;
} Code;

static void build_codes_dfs(Node *n, Code table[256], uint32_t path, uint8_t depth) {
    if (!n) return;
    if (n->is_leaf) {
        table[n->symbol].code = (depth ? path : 0);  /* if single symbol, depth might be 0 */
        table[n->symbol].len  = (depth ? depth : 1); /* ensure at least 1-bit code */
        return;
    }
    /* Left = 0, Right = 1 */
    build_codes_dfs(n->left,  table, (path << 1) | 0u, depth + 1);
    build_codes_dfs(n->right, table, (path << 1) | 1u, depth + 1);
}

static void build_code_table(Node *root, Code table[256]) {
    for (int i = 0; i < 256; ++i) {
        table[i].code = 0;
        table[i].len = 0;
    }
    if (root) build_codes_dfs(root, table, 0, 0);
}

/* -------------------- File Header IO -------------------- */

static const uint8_t MAGIC[4] = { 'H', 'U', 'F', '1' };

static void write_header(FILE *out, uint64_t original_size, const uint32_t freq[256]) {
    fwrite_or_die(MAGIC, 1, 4, out, "fwrite(magic)");
    fwrite_or_die(&original_size, sizeof(original_size), 1, out, "fwrite(size)");
    for (int i = 0; i < 256; ++i) {
        fwrite_or_die(&freq[i], sizeof(uint32_t), 1, out, "fwrite(freq)");
    }
}

static void read_header(FILE *in, uint64_t *original_size, uint32_t freq[256]) {
    uint8_t magic[4];
    if (fread(magic, 1, 4, in) != 4) die_msg("Invalid or truncated file (magic).");
    if (memcmp(magic, MAGIC, 4) != 0) die_msg("Not a HUF1 file (bad magic).");
    if (fread(original_size, sizeof(*original_size), 1, in) != 1) die_msg("Truncated header (size).");
    for (int i = 0; i < 256; ++i) {
        if (fread(&freq[i], sizeof(uint32_t), 1, in) != 1) die_msg("Truncated header (freq).");
    }
}

/* -------------------- Compression -------------------- */

static void compress_file(const char *inpath, const char *outpath) {
    FILE *in = fopen(inpath, "rb");
    if (!in) die("fopen input");
    FILE *out = fopen(outpath, "wb");
    if (!out) die("fopen output");

    uint32_t freq[256] = {0};
    uint64_t original_size = 0;

    /* 1) Count frequencies */
    {
        uint8_t buf[1<<15];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), in)) > 0) {
            original_size += n;
            for (size_t i = 0; i < n; ++i) {
                freq[buf[i]]++;
            }
        }
        if (ferror(in)) die("fread input");
    }

    /* 2) Build tree and code table */
    Node *root = build_huffman_tree(freq);
    Code table[256];
    build_code_table(root, table);

    /* 3) Write header */
    write_header(out, original_size, freq);

    /* 4) Encode data */
    BitWriter bw;
    bw_init(&bw, out);
    if (original_size > 0) {
        if (fseek(in, 0, SEEK_SET) != 0) die("fseek");
        uint8_t buf[1<<15];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), in)) > 0) {
            for (size_t i = 0; i < n; ++i) {
                Code c = table[buf[i]];
                /* Safety: ensure there's a code (length > 0) */
                if (c.len == 0) {
                    /* happens only if freq was zero, which shouldn't occur for seen bytes */
                    die_msg("Internal error: missing code.");
                }
                bw_write_bits(&bw, c.code, c.len);
            }
        }
        if (ferror(in)) die("fread input");
    }
    bw_flush(&bw);

    free_tree(root);
    fclose(in);
    fclose(out);
}

/* -------------------- Decompression -------------------- */

static void decompress_file(const char *inpath, const char *outpath) {
    FILE *in = fopen(inpath, "rb");
    if (!in) die("fopen input");
    FILE *out = fopen(outpath, "wb");
    if (!out) die("fopen output");

    uint64_t original_size = 0;
    uint32_t freq[256];
    read_header(in, &original_size, freq);

    Node *root = build_huffman_tree(freq);

    /* Special cases */
    if (original_size == 0) {
        /* Nothing to decode */
        free_tree(root);
        fclose(in);
        fclose(out);
        return;
    }
    if (!root) {
        die_msg("Corrupt file: missing Huffman tree.");
    }

    BitReader br;
    br_init(&br, in);

    /* If only one distinct symbol existed, tree may have only left child.
       We handle uniformly by walking until we hit a leaf. */
    uint64_t written = 0;
    Node *cur = root;
    int ok = 1;

    while (written < original_size) {
        if (cur->is_leaf) {
            uint8_t byte = (uint8_t)cur->symbol;
            fwrite_or_die(&byte, 1, 1, out, "fwrite(decode)");
            written++;
            cur = root;
            continue;
        }
        int bit = br_read_bit(&br, &ok);
        if (!ok) {
            /* Ran out of bits before writing expected bytes */
            die_msg("Unexpected end of encoded data.");
        }
        cur = (bit == 0) ? cur->left : cur->right;
        if (!cur) die_msg("Corrupt Huffman tree or data.");
    }

    free_tree(root);
    fclose(in);
    fclose(out);
}

/* -------------------- CLI -------------------- */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s -c <input> <output.huf>   Compress\n"
        "  %s -d <input.huf> <output>   Decompress\n",
        prog, prog);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "-c") == 0) {
        compress_file(argv[2], argv[3]);
    } else if (strcmp(argv[1], "-d") == 0) {
        decompress_file(argv[2], argv[3]);
    } else {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}