<?php

/**
 * SudokuPuzzle class.
 *
 * Sudoku puzzle solver that takes in an array (rows) of arrays (columns) of int
 * values (cells) and solves it as a Sudoku puzzle by iterating through a series
 * of checks, each one iterating through each cell in the grid and tracking all
 * possible values for each one, and eventually (hopefully) returning an array
 * of the same structure as the original the represents the puzzle solution.
 *
 *
 * # SOLVING ALGORITHM #
 *
 * We keep track of the puzzle in three variables: `puzzleOrig` (original grid),
 * (possibilities for each cell in the grid -- an array of array of arrays), and
 * `puzzleSol` (solved grid -- or the solution so far).
 *
 * First, for each empty cell, map all possible values (candidates) by looking
 * at each of the other cells in the same row/column/region (each context will
 * be covered individually) and, if our cell's index is another cell's value, we
 * eliminate that candidate from our cell's list.
 *
 * Starting with an array of integers (1-9), cancel out the numbers that appear
 * in the same row (easiest to check), column (values of cells in all other rows
 * in same position as current cell's), and region (starting at the first cell
 * of the region, get each of the three cells in each of the three rows starting
 * with the first cell) as the cell we're iterating through.
 *
 * If this array has a single value, there can be only one possibility, so clear
 * the list of candidates for that cell, and set it to this value. If not, then
 * continue with the next iteration of this loop, if any.
 *
 * After this first sweep, we want to see if any number from 1 to 9 (K) is only
 * possible in ONE cell in any row, column, or region of the puzzle (the "rule
 * of singles"). Basically, we count how many of each of the numbers 1-9 are in
 * the list of possibilities for each empty cell in each area. If a candidate K
 * is is possible in only one of the cells, that cell is the only one that can
 * be that value, so it is then set to that value, and its possibilities list
 * cleared.
 *
 * (TODO: Add some additional algorithms. Keeping it simple [stupid] for now).
 *
 * Finally, we check to see if the puzzle itself is solved. If it is, we return
 * true, ending any calculations and resulting in a solved puzzle grid (which
 * hopefully has NO empty/zero values). If not, we start the series of checks
 * over again (from the top) until either one of two things happen: either the
 * puzzle is eventually solved, at which point we report success; or, it becomes
 * clear that logical deduction will not solve this puzzle. At this point, we
 * stop our looping entirely and try a recursive brute-force algorithm.
 *
 * The brute force algorithm basically goes through each empty cell in the grid,
 * and tries each of its possible candidates, setting the cell to that value as
 * a trial. If this is a valid cell assignment, then we recursively move on to
 * the next empty cell. We repeat this until the puzzle is solved, or we reach
 * the very end of the grid without any success. Each time a placement turns out
 * to be invalid, we report failure (return false), thus backing out of the re-
 * cursion loop one step, which continues the loop and tries the next value, or
 * the next cell if no values are left to try (which is not good).
 *
 */

class SudokuPuzzle {

  ###############
  ## CONSTANTS ##
  ###############

  /**
   * The row index of the puzzle matrix.
   *
   * @var int
   */
  const ROWS = 0;

  /**
   * The column index of the puzzle matrix.
   *
   * @var int
   */
  const COLS = 1;

  /**
   * The region index of the puzzle matrix.
   *
   * @var int
   */
  const BOXS = 2;


  ###############
  ## VARIABLES ##
  ###############

  /**
   * Original Sudoku puzzle, a flat array of 81 (9x9) cells, each represented by
   * an int value.
   *
   * @var array
   */
  private $puzzleOrig = [];

  /**
   * Array of possibilities for each cell index (0-80), each an array of int
   * numbers.
   *
   * @var array
   */
  private $puzzlePoss = [];

  /**
   * Solution to the puzzle so far (or finally), another flat array of ints.
   *
   * @var array
   */
  private $puzzleSol = [];

  /**
   * Number of cells in the grid.
   *
   * @var int
   */
  private $cellCount = 81;

  /**
   * Grid size, number of cells per row/column.
   *
   * @var int
   */
  private $gridSize   = 9;

  /**
   * Region size, number of cells per region, as well as number of regions per
   * row/column overall.
   *
   * @var int
   */
  private $regionSize = 3;

  /**
   * An array of arrays of arrays of cell indexes (integers). Basically, a list
   * of three grids (rows, columns, regions), each of which is an array of rows,
   * each of which is a list of cell IDs in that row, column, or region.
   *
   * @var array
   */
  private $puzzleMatrix = [];


  #################
  ## CONSTRUCTOR ##
  #################

  /**
   * Instantiate and initialize our Sudoku puzzle and solution grids, as well as
   * the initial state of our grid of "Possibilities" arrays.
   *
   * @param array $sudokuGrid
   */
  public function __construct($sudokuGrid = []) {
    // Simple validation, checking for an array of NxN values.
    if (!is_array($sudokuGrid)) {
      throw new Error("Sudoku grid must be an array");
    }
    $gridCount = count($sudokuGrid);
    $gridSize  = sqrt($gridCount);
    $blockSize = sqrt($gridSize);
    if (floor($gridSize) !== $gridSize || floor($blockSize) !== $blockSize) {
      // Given array length not perfect square, or cannot be split into regions.
      throw new Error(
        "Sudoku grid must be array of NxN numbers, where N is a square number."
      );
    }

    // Store our size values.
    $this->cellCount  = $gridCount;
    $this->gridSize   = $gridSize;
    $this->regionSize = $blockSize;

    // Initialize our puzzle matrix, just a helpful list of grids.
    $this->puzzleMatrix = $this->getPuzzleMatrix();

    // Sanitize array by ensuring each is an int between 0 and 9, defaulting
    //  to 0 if false/null/undefined/anything else funky and unexpected.
    $sudokuGrid = array_map(fn($val) => intval($val), $sudokuGrid);

    // Initialize original and solution grids with original puzzle passed in.
    $this->puzzleOrig = $sudokuGrid; // to preserve original puzzle
    $this->puzzleSol  = $sudokuGrid; // iteratively, track solution

    // Next, set possibilities of each EMPTY cell to any value - i.e. [1,...,9].
    // Each filled cell gets empty array, which we'll skip each solve iteration.
    for ($N = 0; $N < $this->cellCount; $N++) {
      $this->puzzlePoss[$N] = ($this->puzzleSol[$N] > 0)
        ? [] : $this->puzzlePoss[$N] = range(1, 9);
    }
  }


  #########################################
  ## PUBLIC METHODS / MAIN FUNCTIONALITY ##
  #########################################

  /**
   * Called when we want to solve a given puzzle.
   *
   * @return bool
   */
  public function solve(): bool {
    // Dangerous, but easier. We should break free eventually no matter what.
    while (true) {

      // Keep track of whether we were able to make ANY changes on this pass.
      $previousState = $this->puzzleSol;

      // First, [re-]calculate the possibilities for each cell in solution grid.
      $this->processPossibsForEachCell();

      // Next, check for any cells that can be solved by the rule of singles.
      $this->processRuleOfSingles();

      // Next, check for any groupings of cells whose candidates are exclusive.
      $this->processNakedSubsetRule();

      // After all of the above, do a final check for single-candidate arrays.
      $this->processSingleCandidates();

      // If we've solved the puzzle already, there's nothing more to do.
      if ($this->isSolved()) return true;

      // If no changes were made in this iteration, exit loop for another pass.
      if ($previousState === $this->puzzleSol) break;

    } # while true

    // Implement brute force recursively. We've exhausted all logical means.
    return $this->recursiveBruteForce(); // by default starts with cell 0
  }

  /**
   * Algorithm Passthrough, Step 1
   *
   * Determine, or recalculate, the list of possible values for each cell in the
   * Sudoku solution grid. Run on each iterative pass of the `solve` algorithm.
   *
   * @return void
   */
  private function processPossibsForEachCell(): void {
    for ($N = 0; $N < $this->cellCount; $N++) {
      $this->puzzlePoss[$N] = array_values($this->puzzlePoss[$N]);
      // Run check on possibilities for each cell N in case any must be removed.
      $this->calcPossibsForCell($N);
    }

    // Should be done at the end of every passthrough step function.
    $this->performCandidateSweep();
  }

  /**
   * Algorithm Passthrough, Step 2
   *
   * Here, we want to see if any of the candidate values from 1 to 9 (K) are
   * possible in only ONE empty cell of each row, column, and 3x3 region in the
   * puzzle solution  (i.e. "rule of singles"). If so, that cell MUST contain
   * that value K. Basically, loop through the cells in each row, column, and
   * region, then, for each cell N in each of the 3 lists, start a counter to
   * keep track of how many of cells in that list contain in their respective
   * lists of possible candidates any of the values K, and how many there are.
   * Iterate through each of number 1 through 9: for each number K, reset the
   * counter and increment it each time we find a cell whose list of possible
   * candidates (`puzzlePoss[N]`) contains K. If this counter is equal to ONE,
   * we set the value of cell N to this number K, then clear the array of that
   * cell's possibilites entirely (`puzzlePoss[N] = []`).
   *
   * @return void
   */
  private function processRuleOfSingles(): void {
    $matrix = $this->puzzleMatrix; // stored in constructor

    // Check the list of IDs for each row (index 0), column (1), and region (2).
    foreach ($matrix as $gridType => $grid) {
      $gridTypeStr = $gridType < 1 ? 'ROW' : ($gridType < 2 ? 'COL' : 'REG');

      // Check the list of IDs for this row/column/region.
      foreach ($grid as $listIdx => $list) {

        // Check each possible candidate value K from 1 to 9.
        foreach (range(1, $this->gridSize) as $dblK) {
          $K = intval($dblK); // not sure why these are doubles, but OH WELL
          // If this candidate is already known in this row/column/region, skip.
          foreach ($list as $checkCellId) {
            if ($this->puzzleSol[$checkCellId] === $K) continue 2;
          }

          // Initialize counter used to track # of cells with K as a candidate.
          $counter = 0;

          // Check each cell ID in the list for this particlar row/column/box.
          foreach ($list as $cellId) {
            // Check the value of the cell in this ID (0 means unknown).
            $val = intval($this->puzzleSol[$cellId]);

            // Skip if this already known.
            if ($val > 0) continue;

            $this->puzzlePoss[$cellId] = array_map('intval', $this->puzzlePoss[$cellId]);
            // If this cell is unknown and K is one of its candidates, count it.
            if ($val === 0 && in_array($K, $this->puzzlePoss[$cellId], true)) {
              // Store cell ID in case there's only 1 cell with K as candidate.
              $lastCellId = $cellId;
              // Mark this as one [more] cell that could possibly contain K.
              $counter++;
            }

          } # each $cellId

          // If K is only possible is ONE cell of this grouping, solve for cell.
          if ($counter === 1) {
            // If one and only one, we've solved for this sole cell (equals K).
            $this->solveCell($lastCellId, $K);
          }

        } # each $K

      } # each $list

    } # each $grid

    // Should be done at the end of every passthrough step function.
    $this->performCandidateSweep();

  }

  /**
   * Algorithm Passthrough, Step 3
   *
   * The Naked Subset Rule states that if N candidates are possible in a certain
   * set of N cells, all in the same row, column, or region, and no other candi-
   * dates are possible in those cells, then those N candidates are not possible
   * elsewhere in that same row, column, or region.
   *
   * For example, if two empty cells in a particular row/column/region had 3 and
   * 5 as the only two possibilities for both, then all of the other empty cells
   * must have 3 and 5 REMOVED from their list of possible values, as only those
   * 2 cells can have either of those 2 numbers. The same would apply to a group
   * of 3 cells each with the same list of 3 candidates, and so on.
   *
   * @return void
   */
  private function processNakedSubsetRule(): void {
    // Store shortcuts.
    $solution = $this->puzzleSol;
    $matrix   = $this->puzzleMatrix;

    // Go through each set of rows, columns, and regions.
    // Each grid is an array (row/column/block) of arrays (cell IDs to check).
    foreach ($matrix as $gridType => $grid) {
      $gridTypeName = $gridType === self::ROWS ? 'row'
        : ($gridType === self::COLS ? 'column' : 'region');
      // echo "<h1>Going through each {$gridTypeName}...</h1>";

      // "Area" in this case means each row, column, or region in the loop.
      foreach ($grid as $areaId => $area) {
        // echo "<h2>", ucfirst($gridTypeName), " #{$areaId}...</h2>";

        // Each row/column/region entry is a list of all 9 cells in that area.
        foreach ($area as $cellId) {

          // This cell we're iterating through will be referred to as cell "N".
          $cellN = $solution[$cellId];
          // We only want unknown/empty cells.
          if ($cellN > 0) continue;

          // echo "<h3>Empty cell (# {$cellId})...</h3>";

          // Keep track of the cells we've already checked as part of a subset.
          $checked = [];
          // Store all possible values for cell N.
          $possN = $this->puzzlePoss[$cellId];

          // echo "<h4>Possibilites: [", implode(', ', $possN), "]</h4>";

          // Start our counter at 1 to include cell N as first.
          $counter = 1;
          // Count the current cell as the first one checked.
          $checked[] = $cellId;

          // Go through all of the other cells and check their candidates lists.
          foreach ($area as $otherCellId) {

            // Skip the current cell (we want all *other* cells).
            if ($cellId === $otherCellId) continue;

            // This cell will be known as cell "M"; grab its candidates list.
            $possM = $this->puzzlePoss[$otherCellId];

            // echo "<h5>Checking other cell in {$gridTypeName} → ",
            //      "index {$otherCellId}...</h5>";

            // Check if cell M's candidates list is identical to cell N's list.
            if (array_values($possN) === array_values($possM)) {
              // Count this as part of the subset (so far).
              $counter++;

              // echo "<p style='color:blue'>Found an exact match for candidates",
              //      " list! Number of cells so far = {$counter}.</p>";

              // Count it as already checked, to avoid duplicate checks.
              $checked[] = $otherCellId;
            }

          } # each $otherCellId

          // If our subset doesn't count as a naked one, skip it and move on.
          $possNcount = count($possN);
          if ($counter !== $possNcount) {
            continue;
          }

          // Otherwise, naked subset -- no longer possible for any of the given
          // cells to have as their value any candidate in our subset.
          foreach ($area as $otherId) {
            // This isn't an "other cell", so skip it.
            if (in_array($otherId, $checked)) continue;

            // Go through each candidate in the cells' shared list.
            foreach ($possN as $K) {
              // If this cell doesn't have this candidate in its list, skip.
              if (!in_array($K, $this->puzzlePoss[$otherId])) continue;

              // Candidate K can't be the value of this cell outside our subset.
              $idxOfK  = array_search($K, $this->puzzlePoss[$otherId]);
              $cellVal = $this->puzzlePoss[$otherId];
              $this->puzzlePoss[$otherId] = array_splice($cellVal, $idxOfK, 1);
            } # each $K

          } # each $otherId

        } # each $cellId

      } # each $area

    } # each $grid

    // This should be done at the end of every passthrough step function.
    $this->performCandidateSweep();

  }

  /**
   * Algorithm Passthrough, Step 4
   *
   * After all of the above, we want to check for any lists of candidates that
   * end up being a single value. If this is the case, we mark that cell as
   * solved for that single value and empty out the candidate possibilities
   * list.
   *
   * @return void
   */
  private function processSingleCandidates(): void {
    // For each possibilities array, if there is only one candidate, solve cell.
    foreach ($this->puzzlePoss as $N => $possN) {
      // Check for single-candidate array.
      if (count($possN) === 1) {
        // Mark cell N as solved, with single value stored in candidates array.
        $this->solveCell($N, $possN[0]);
      }
    }
  }

  /**
   * Algorithm Passthrough -- End of Each Step
   *
   * After each step above, we want to check for any lists of candidates that
   * end up being a single value. If this is the case, we mark that cell as
   * solved for that single value and empty out the candidate possibilities
   * list.
   *
   * @return void
   */
  private function performCandidateSweep(): void {
    $candidates = $this->puzzlePoss;
    // Go through the candidates for each cell in the solution grid.
    foreach ($candidates as $cellId => $list) {
      // If there is only 1 possibility, logically that cell must be that value.
      if (count($list) === 1) {
        // Set this cell's value (also clears its list of possible candidates).
        $this->solveCell($cellId, $list[0]);
      }
    }
  }

  /**
   * Fallback for Deduction-Based Algorithms -- Recursive Brute Force
   *
   * Brute force recursive algorithm. Basically, we go from empty cell to empty
   * cell and try each candidate in its list. Then, we keep going, doing the
   * same thing until we either: solve the puzzle (should happen eventually); or
   * fail, and thus fall back to the next candidate in the list of the cell that
   * failed. This loop cycles repeatedly until the puzzle has been solved by
   * brute force. Used when all other methods have failed to solve the puzzle.
   *
   * @return bool Success?
   */
  private function recursiveBruteForce($cellId = 0): bool {
    // Base case: Reached past end of Sudoku puzzle grid.
    if ($cellId >= $this->cellCount) return true;
    // If this cell is already solved, skip over it to the next one.
    if ($this->puzzleSol[$cellId] > 0) {
      return $this->recursiveBruteForce($cellId + 1);
    }
    // Try each number in this cell's list of candidates.
    $tryList = $this->puzzlePoss[$cellId]; // separate, so we can alter as we go
    foreach ($tryList as $K) {
      if ($this->isSafeToPlace($cellId, $K)) {
        // Try this value on for size.
        $this->puzzleSol[$cellId] = $K;
        // Check to see if this works.
        if ($this->recursiveBruteForce($cellId + 1)) {
          return true;
        }
        // Nope, that didn't work. Back out and try the next candidate, or cell.
        $this->puzzleSol[$cellId] = 0;
      }
    }
    // Nothing worked. Report utter failure.
    return false;
  }

  /**
   * Check the candidates for the cell at index N against the list of cells
   * passed in. If any of the cells contain the candidate K, we remove K from
   * the list of candidates for cell N. If K is found in any of the cells, we
   * return true, indicating that K was culled from the list of candidates.
   *
   * @param integer $K Candidate value to check against cells in `$cells` list.
   * @param array $cells Array of all cells in the row/column/region.
   * @param integer $N Index of cell whose candidate list to check against.
   * @return boolean Did we end up having to cull K from the list of candidates?
   */
  private function cullCandidate(int $K, array $cells, int $N): bool {
    // Basically, we need to check the cells represented by each of the indexes
    // in array `cells` and, if any contain the given candidate K, we remove the
    // candidate K from the list of possibilities (at K index idxK) for the cell
    // at index N.
    foreach ($cells as $cellId) {
      if ($cellId === $N) continue; // skip self
      $cellVal = $this->puzzleSol[$cellId];
      if ($cellVal === $K) {
        return true; // no need to check any more cells (we found K)
      }
    }
    // If we made it to the end, we must not have found it, so we didn't cull.
    return false; // this, strangely, means "success" vàv remaining a candidate
  }


  ############################
  ## GRID UTILITY FUNCTIONS ##
  ############################

  /**
   * If marking the cell at `$cellId` with value `$value` doesn't violate any
   * rules, then it is deemed safe to place, at least initially.
   *
   * @param integer $cellId Cell index (0-80).
   * @param integer $value Value to attempt to place.
   * @return boolean Is it safe to place this value here?
   */
  private function isSafeToPlace(int $cellId, int $value): bool {
    $gridSize = $this->gridSize;
    $regSize  = $this->regionSize;

    // Calculate the row and column indices based on the cell ID.
    $thisRowId = floor($cellId / $gridSize);
    $thisColId = $cellId % $gridSize;

    // Test all cells in the same row.
    for ($testCol = 0; $testCol < $gridSize; $testCol++) {
      if ($this->puzzleSol[$thisRowId * $gridSize + $testCol] === $value) {
        return false;
      }
    }

    // Test all cells in the same column next.
    for ($testRow = 0; $testRow < $gridSize; $testRow++) {
      if ($this->puzzleSol[$testRow * $gridSize + $thisColId] === $value) {
        return false;
      }
    }

    // Test all cells in the same region.
    $startRow = $thisRowId - ($thisRowId % $regSize);
    $startCol = $thisColId - ($thisColId % $regSize);
    for ($rowOff = 0; $rowOff < $regSize; $rowOff++) {
      for ($colOff = 0; $colOff < $regSize; $colOff++) {
        $testCellId = ($startRow + $rowOff) * $gridSize + ($startCol + $colOff);
        if ($this->puzzleSol[$testCellId] === $value) {
          return false;
        }
      }
    }

    // If no conflicts were found, it is safe to place the value.
    return true;
  }

  /**
   * Mark the given cell as having a known value (also given). This clears that
   * cell's possibilities array, as it is no longer unknown. The wave function
   * has collapsed, so to speak.
   *
   * @param integer $cellId Cell index (0-80).
   * @param integer $value Value to place in this cell.
   * @return void
   */
  private function solveCell(int $cellId, int $value): void {
    $this->puzzleSol[$cellId]  = $value;
    $this->puzzlePoss[$cellId] = [];
  }

  /**
   * Check each candidate in the list of possibilities for cell N to see if any
   * of the other cells in the same row, column, or region contain them. Remove
   * any candidate found, thus eliminating it from the list of possible values.
   *
   * @param integer $N Cell index to check.
   * @return void
   */
  private function calcPossibsForCell(int $N): void {
    if (empty($this->puzzlePoss[$N])) return;

    // Gather our cell's row, column, and region indexes (all 0-8).
    $cellRow = $this->getRowIndex($N);
    $cellCol = $this->getColIndex($N);
    $cellReg = $this->getRegionIndex($N);

    // We run our checks on this temporary array, but we alter the original.
    $possTmp = array_values($this->puzzlePoss[$N]);

    // Go through each remaining candidate, and remove if found anywhere.
    $newPossArray = []; // candidates (K) that we're KEEPING
    foreach ($possTmp as $indexOfK => $K) {

      // First, check cells in the same row as the current cell.
      $otherCells = $this->getCellIdsInRow($cellRow);
      // If any of these cells contain K as a value, remove from possibilities.
      $wasFound = $this->cullCandidate($K, $otherCells, $N);
      // If we found K, no need to keep looking.
      if ($wasFound) continue;

      // Next, check cells in the same column as the current cell.
      $otherCells = $this->getCellIdsInCol($cellCol);
      // If any of these cells contain K, remove K from list of possibilities.
      $wasFound = $this->cullCandidate($K, $otherCells, $N);
      // If K found, stop looking.
      if ($wasFound) continue;

      // Finally, check cells in same region.
      $otherCells = $this->getCellIdsInReg($cellReg);
      // If any contain K, eliminate K as candidate.
      $wasFound = $this->cullCandidate($K, $otherCells, $N);
      if ($wasFound) continue;

      // If we've gotten this far, K is still a candidate; add to the keep list.
      $newPossArray[] = $K;
    }

    // If the new array has only ONE element, this cell must contain that value.
    if (count($newPossArray) === 1) {
      $this->solveCell($N, $newPossArray[0]);
      return;
    }

    // [Re-]Set the possibilities for this cell.
    $this->puzzlePoss[$N] = $newPossArray;
  }

  /**
   * Reset the array keys of each empty cell's list of candidates (done after we
   * cull candidates from any of the cells' lists).
   *
   * @return void
   */
  private function resetPossibilitiesArrayKeys(): void {
    // We look at each cell's list of candidates (known cells will be skipped).
    foreach ($this->puzzlePoss as $N => $cellPoss) {
      // We need to reset the keys of this empty cell's candidates.
      $this->puzzlePoss[$N] = array_values($cellPoss);
    }
  }


  ############################
  ## DATA UTILITY FUNCTIONS ##
  ############################

  /**
   * Return the list of cell indices within the specified row/column/region (the
   * type of which is also passed in).
   *
   * @param integer $areaId Index of the row, column, or region to look at.
   * @param integer $areaType Type of area (i.e. if a row, column, or region).
   * @return array List of each index of all cells in the given row/col./region.
   */
  private function getCellIdsInArea(int $areaId, int $areaType): array {
    switch ($areaType) {
      // If the given index is of a row, get all cells in that row.
      case self::ROWS:
        return $this->getCellIdsInRow($areaId);
      // If the given index is of a column, get all cells in that column.
      case self::COLS;
        return $this->getCellIdsInCol($areaId);
      // If the given index is of a region, get all cells in that region.
      case self::BOXS;
        return $this->getCellIdsInReg($areaId);
      // If anything else, we don't know what to check so NO cells (empty list).
      default:
        return [];
    }
  }

  /**
   * Return the list of cell indices in the given row.
   *
   * @param integer $rowIdx Index of the row to look at.
   * @return array List of cells in the given row.
   */
  private function getCellIdsInRow(int $rowIdx): array {
    $firstCellIdx = $rowIdx * $this->gridSize;
    // Generate array of ints starting at first cell in row, step by 1, 9 times.
    return $this->range($firstCellIdx, 1, $this->gridSize);
  }

  /**
   * Return the list of cell indices in the given column.
   *
   * @param integer $colIdx Index of the column to look at.
   * @return array List of cells in the given column.
   */
  private function getCellIdsInCol(int $colIdx): array {
    // Generate array of ints starting at col N (which is the same as the Nth
    //  cell in the first row), incrementing by 9 (i.e. jumping rows), 9 times.
    return $this->range($colIdx, $this->gridSize, $this->gridSize);
  }

  /**
   * Return the list of cell indices in the given region.
   *
   * @param integer $regIdx Index of the region to look at.
   * @return array List of cells in the given column.
   */
  private function getCellIdsInReg(int $regIdx): array {
    // Generate array of ints that lie within region N (0-8 in a 3-by-3 grid).
    if ($regIdx < 0 || $regIdx > $this->regionSize ** 2) return [];
    // Start with the first cell in the region.
    $firstId = (intdiv($regIdx, 3) * 27) + (($regIdx % 3) * 3);
    $cellIds = []; // empty - we'll fill in first cell along with the rest
    $offsets = [0, 1, 2];
    // Go through each of the three rows in this region.
    foreach ($offsets as $rowOff) {
      // Go through each of the three columns in this region.
      foreach ($offsets as $colOff) {
        // Calculate this particular cell index based on our offset values.
        $thisId = $firstId + (($this->regionSize ** 2) * $rowOff) + $colOff;
        // Store this cell index.
        $cellIds[] = $thisId;
      }
    }
    // Since we've already dealt with the skip condition (if any) return result.
    return $cellIds;
  }

  /**
   * Return the index of the row in which the given cell lies.
   *
   * @param integer $N Cell index to look for.
   * @return integer Row index of this cell.
   */
  private function getRowIndex(int $N): int {
    if (!is_numeric($N) || !$this->isValidCellIdx($N)) return false;
    return floor($N / $this->gridSize); // e.g. cell idx 47 / 9 ==> row idx 5
  }

  /**
   * Return the index of the column in which the given cell lies.
   *
   * @param integer $N Cell index to look for.
   * @return integer Column index of this cell.
   */
  private function getColIndex(int $N): int {
    if (!is_numeric($N) || !$this->isValidCellIdx($N)) return false;
    return $N % $this->gridSize; // e.g. cell idx 47 % 9 ==> col idx 2
  }

  /**
   * Return the index of the region in which the given cell lies.
   *
   * @param integer $N Cell index to look for.
   * @return integer Region index of this cell.
   */
  private function getRegionIndex(int $N): int {
    if (!is_numeric($N) || !$this->isValidCellIdx($N)) return false;
    $regnSize = $this->regionSize;
    $rowN = $this->getRowIndex($N);
    $colN = $this->getColIndex($N);
    if ($rowN === false || $colN === false) return false;
    /**
     * RegionIndex = 3 x   RegionRowIndex   +   RegionColIndex
     *             = 3 x ⌊CellNRowIndex / 3⌋ + ⌊CellNColIndex / 3⌋
     */
    return $regnSize * intdiv($rowN, $regnSize) + intdiv($colN, $regnSize);
  }

  /**
   * Return whether the given cell index is a valid one for this grid (i.e. it
   * is between 0 and the highest cell index possible).
   *
   * @param integer $N Cell index to check.
   * @return boolean Is this a valid cell index?
   */
  private function isValidCellIdx(int $N): bool {
    // N must be between 0 and the highest cell index possible, inclusive.
    return $N >= 0 && $N < $this->cellCount;
  }


  ##############################
  ## OUTPUT UTILITY FUNCTIONS ##
  ##############################

  /**
   * Check if the puzzle is solved, where every cell's value is greater than 0.
   *
   * @param boolean $asText Whether to output text (or array).
   * @return array|string
   */
  public function getPuzzle($asText = true): array|string {
    return $asText ? implode('', $this->puzzleOrig) : $this->puzzleOrig;
  }

  /**
   * Return the solution to the puzzle.
   *
   * @param boolean $asText Whether to output text (or array).
   * @return array|string
   */
  public function getSolution($asText = true): array|string {
    return $asText ? implode('', $this->puzzleSol) : $this->puzzleSol;
  }

  /**
   * Return the puzzle's solution grid as HTML, for debugging purposes.
   *
   * @return string Solution grid HTML.
   */
  public function getSolutionHtml(): string {
    $html = '<div class="sudoku">';
    // Build by region; makes the most sense HTML-wise.
    for ($region = 0; $region < $this->regionSize ** 2; $region++) {
      $html .= '<div class="region">';
      $cellsInRegion = $this->getCellIdsInReg($region);
      foreach ($cellsInRegion as $cellId) {
        $sol  = $this->puzzleSol[$cellId] > 0 ? $this->puzzleSol[$cellId] : '';
        $orig = ($this->puzzleOrig[$cellId] === $this->puzzleSol[$cellId]);
        $html .= '<div class="cell' . ($orig ? '' : ' new') . "\">{$sol}</div>";
      }
      $html .= '</div>';
    }
    $html .= '</div>';
    return $html;
  }

  /**
   * Return the SVG representation of the Sudoku puzzle as an XML string.
   *
   * @param boolean $getSolution Should we get the solution? If not, get orig.
   * @param boolean $showSolvedDigits Show solved digits as darker? N/A if orig.
   * @return string Sudoku grid SVG/XML.
   */
  public function getPuzzleSvg(
    $getSolution = false,
    $showSolvedDigits = false
  ): string {
    // Load the SVG template.
    $empty_svg = new DOMDocument();
    $empty_svg->load('images/sudoku/sudoku_grid.svg');
    // Get root element of SVG document.
    $svg = $empty_svg->documentElement;
    // Iterate over each text#Digit_N element in the SVG source.
    $elements = $svg->getElementsByTagName('text');
    foreach ($elements as $index => $element) {
      // Check ID and grab corresponding value from puzzle/solution.
      $element_id = $element->getAttribute('id');
      // Get the cell index from the element ID.
      $cell_id = str_replace('Digit_', '', $element_id);
      // Get the digit from the
      $digit = $getSolution === true // only if explicitly true
        ? $this->puzzleSol[$cell_id] : $this->puzzleOrig[$cell_id];
      $digit = $digit > 0 ? "$digit" : ''; // show unknown as blank
      // Set text content of this element to the digit found.
      $element->textContent = $digit;
      // If this is a cell that we solved, and we're told to, show as green.
      if (
        $showSolvedDigits === true && $getSolution === true &&
        $this->puzzleSol[$cell_id] !== $this->puzzleOrig[$cell_id]
      ) {
        $element->setAttribute('class', 'new');
      }
    }
    // Minimize the SVG.
    $empty_svg->formatOutput = false;
    $empty_svg->preserveWhiteSpace = false;
    $empty_svg->normalizeDocument();
    // Output the modified SVG.
    return $empty_svg->saveXML();
  }

  /**
   * Return a three-tiered matrix consisting of a list of 3 grids: list of cell
   * IDs in each row, then each column, then each 3x3 region.
   *
   * @return array Matrix of cells to check each pass, organized by type/region.
   */
  private function getPuzzleMatrix(): array {
    $size = $this->gridSize;
    $matrix = [ [], [], [] ];
    // Index 0 - Rows.
    for ($i = 0; $i < $size; $i++) {
      // The Ith row of the ROWS grid is a list of cell IDs in that row.
      $matrix[self::ROWS][$i] = $this->getCellIdsInRow($i);
    }
    // Index 1 - Columns.
    for ($i = 0; $i < $size; $i++) {
      // The Ith column of the COLS grid is a list of cell IDs in that column.
      $matrix[self::COLS][$i] = $this->getCellIdsInCol($i);
    }
    // Index 2 - Regions/boxes.
    for ($i = 0; $i < $size; $i++) {
      // The Ith region of the BOXS grid is a list of cell IDs in that region.
      $matrix[self::BOXS][$i] = $this->getCellIdsInReg($i);
    }
    // Aaaand we're done.
    return $matrix;
  }

  /**
   * Whether the puzzle is currently in a solved state (i.e. no 0 values left).
   *
   * @return bool Whether puzzle has been / is solved
   */
  private function isSolved(): bool {
    // Return whether the list of all 0s in the solution is empty (i.e. solved).
    return empty(array_filter($this->puzzleSol, fn($val) => $val === 0));
  }

  /**
   * Semi-alias for PHP's `range` that generates an array of integers starting
   * at a given value, incrementing by a specified value, and of a given length.
   *
   * @param int $start Starting value of the sequence.
   * @param int $step Increment value between each integer.
   * @param int $length Total number of integers in the array.
   * @return array Array of integers.
   */
  private function range(int $start = 0, int $step = 1, int $len = 0): array {
    if ($len <= 0) return [];
    return range($start, $start + ($len - 1) * $step, $step);
  }

  /**
   * Takes in an array and removes one or all element(s) matching `$value`, opt-
   * ionally preserving the original keys.
   *
   * @param array $source Original array
   * @param mixed $value Value to remove in original array
   * @param bool $saveKeys Preserve the keys of the original array?
   * @param bool $removeOne Remove only one element (true) or all (false)?
   * @return array
   */
  private function arrayWithout(
    array $source    = [],
    mixed $value     = null,
    bool  $saveKeys  = true,
    bool  $removeOne = false
  ): array {
    // Check for validity of both source and value to remove.
    if (empty($source) || !is_array($source)) return []; // [empty] - ? = empty
    if ($value === null || !isset($value)) return [];  // ? - [invalid] = empty

    // Keep tracking of our working array in a new variable, just in case.
    $newArray = $source;

    // Get an array of keys to remove, even if only one.
    $keysToRemove = ($removeOne === true)
      ? [array_search($value, $newArray)] // returns first key, so array of one
      :  array_keys($newArray, $value); // returns ALL keys where elem == value

    // Remove all keys (or the one key) we need to remove our working array.
    foreach ($keysToRemove as $key) {
      unset($newArray[$key]);
    }

    // If we're explicitly not preserving keys, just return a flat array.
    return ($saveKeys === false)
      ? array_values($newArray) // flattens array
      : $newArray;  // preserves keys from source
  }

}
