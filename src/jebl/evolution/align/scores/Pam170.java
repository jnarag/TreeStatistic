package jebl.evolution.align.scores;

public class Pam170 extends AminoAcidScores {

  private final float[][] residueScores = {

            /*  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V */
            {   3},
            {  -3,  8},
            {   0, -1,  4},
            {   0, -3,  3,  6},
            {  -3, -5, -5, -7, 13},
            {  -1,  1,  0,  2, -8,  6},
            {   0, -2,  2,  5, -8,  3,  6},
            {   1, -4,  0,  0, -5, -2,  0,  6},
            {  -3,  2,  2,  0, -5,  4,  0, -4,  9},
            {  -1, -3, -3, -4, -3, -3, -3, -4, -4,  7},
            {  -3, -4, -4, -6, -9, -2, -5, -6, -3,  2,  7},
            {  -2,  4,  1, -1, -8,  0, -1, -3, -1, -3, -4,  6},
            {  -2, -1, -3, -4, -7, -1, -3, -4, -4,  2,  4,  1, 10},
            {  -5, -6, -5, -8, -6, -7, -8, -6, -3,  1,  1, -8,  0, 10},
            {   1, -1, -1, -2, -4,  0, -1, -2, -1, -3, -4, -2, -3, -6,  8},
            {   2, -1,  1,  0,  0, -1, -1,  1, -2, -2, -4, -1, -2, -4,  1,  3},
            {   2, -2,  0, -1, -3, -2, -1, -1, -2,  0, -3,  0, -1, -5,  0,  2,  5},
            {  -8,  2, -5, -9,-10, -7,-10, -9, -4, -7, -3, -5, -6, -1, -8, -3, -7, 18},
            {  -5, -6, -3, -6,  0, -6, -6, -7,  0, -2, -2, -6, -4,  7, -7, -4, -4, -1, 12},
            {   0, -4, -3, -4, -3, -3, -3, -2, -3,  5,  2, -4,  2, -2, -2, -2,  0, -9, -4,  6}};

  public Pam170() { buildScores(residueScores); }
}
