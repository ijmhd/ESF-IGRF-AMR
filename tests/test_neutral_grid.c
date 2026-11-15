#include "neutral_grid.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "neutral grid failed line %d\n", __LINE__);                            \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
/**
 * Write a small valid neutral-grid fixture with exact header, position binding and checksum.
 */
static int write_grid(const char *path, TubeGrid *g, double v[12]) {
    FILE *f = fopen(path, "wb");
    if (!f)
        return -1;
    uint64_t pos = UINT64_C(14695981039346656037);
    for (int i = 0; i < 2; i++)
        pos = neutral_grid_hash(pos, &g->cell[i].position, 24);
    uint64_t h[] = {1, 2, 1, 1, 6, pos, neutral_grid_hash(UINT64_C(14695981039346656037), v, 96)};
    double t[] = {1000, 1500};
    int ok = fwrite("ESFNEU01", 1, 8, f) == 8 && fwrite(h, 8, 7, f) == 7 &&
             fwrite(t, 8, 2, f) == 2 && fwrite(v, 8, 12, f) == 12;
    return fclose(f) || !ok ? -1 : 0;
}
/**
 * Check spatial neutral format, geometry/temperature/hash binding and corruption rejection.
 */
int main(void) {
    const char *path = "neutral_grid_test.bin", *csv = "neutral_grid_test.csv";
    TubeCell cells[2] = {{.height = 300e3,
                          .position = {6671200, 0, 0},
                          .B = 3e-5,
                          .b = {0, 0, 1},
                          .gx = {1, 0, 0},
                          .gy = {0, 1, 0}},
                         {.height = 350e3,
                          .position = {6721200, 0, 0},
                          .B = 4e-5,
                          .b = {0, 0, 1},
                          .gx = {1, 0, 0},
                          .gy = {0, 1, 0}}};
    TubeCell reference[2];
    memcpy(reference, cells, sizeof(cells));
    TubeGrid g = {.c = {.nx = 2, .ny = 1, .ns = 1}, .cell = cells}, old = g;
    old.cell = reference;
    double v[] = {1e16, 2e15, 3e14, 800, .2, 20, 4e15, 5e14, 6e13, 900, .02, 2}, *data = NULL;
    CHECK(!write_grid(path, &g, v));
    CHECK(neutral_grid_read(path, &g, 1000, 1500, &data) == 1);
    CHECK(!memcmp(data, v, sizeof(v)));
    tube_conductivity_values(&g, 1e10, data);
    free(data);
    FILE *f = fopen(csv, "w");
    CHECK(f);
    fputs("height_m,nu_O_s1,nu_e_s1\n300000,.2,20\n350000,.02,2\n", f);
    CHECK(!fclose(f));
    CHECK(!tube_conductivity_table(&old, 1e10, csv));
    for (int i = 0; i < 2; i++) {
        CHECK(fabs(cells[i].nu / reference[i].nu - 1) < 1e-14);
        CHECK(fabs(cells[i].pF / reference[i].pF - 1) < 1e-14);
        CHECK(fabs(cells[i].hF / reference[i].hF - 1) < 1e-14);
        for (int j = 0; j < 2; j++)
            CHECK(fabs(cells[i].cF[j] - reference[i].cF[j]) < 1e-20);
    }
    CHECK(neutral_grid_read(csv, &g, 1000, 1500, &data) == 0);
    CHECK(neutral_grid_read(path, &g, 999, 1500, &data) == -1);
    cells[0].position.x += 1;
    CHECK(neutral_grid_read(path, &g, 1000, 1500, &data) == -1);
    cells[0].position.x -= 1;
    f = fopen(path, "r+b");
    CHECK(f);
    CHECK(!fseek(f, 87, SEEK_SET));
    int c = fgetc(f);
    CHECK(c != EOF);
    CHECK(!fseek(f, 87, SEEK_SET));
    CHECK(fputc(c ^ 1, f) != EOF);
    CHECK(!fclose(f));
    CHECK(neutral_grid_read(path, &g, 1000, 1500, &data) == -1);
    v[0] = NAN;
    CHECK(!write_grid(path, &g, v));
    CHECK(neutral_grid_read(path, &g, 1000, 1500, &data) == -1);
    remove(path);
    remove(csv);
    puts(
        "PASS spatial response equivalence; mapping, temperatures, corruption and nonfinite rejection");
    return 0;
}
