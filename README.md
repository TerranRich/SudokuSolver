# SudokuSolver
This is **SudokuSolver**, a PHP class that solves nearly any valid 3x3 Sudoku puzzle using both logical deduction (first) and recursive brute force (last resort).

## Usage
After including the `SudokuSolver` class (or the `SudokuSolver.php` file), just pass a flat array of 81 numbers to the constructor. They can be strings, but they will be sanitized to integers upon instantiation.

Then, to solve, use `solve()`, and check the return value. At any point, you can access `getPuzzle()` to retrieve an array representing the original puzzle passed in, and `getSolution` to retrieve an array representing the solution (or the solution so far).

If you need to output the solution as HTML, here are the two variables that I use, with the generated Sudoku grid HTML styled appriopriately:

```php
// If we need HTML (for, say, the partial solution if solving failed), use this as a basis.
$htmlHeader = '<!DOCTYPE html>
  <html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sudoku Solver</title>
  </head>
  <body>
  <style>
    .sudoku {
      border: 2px solid black;
      border-collapse: collapse;
      display: grid;
      font: 18px "Roboto", Arial, Helvetica, sans-serif;
      grid-template-columns: 4.5rem 4.5rem 4.5rem;
      grid-template-rows: 4.5rem 4.5rem 4.5rem;
      user-select: none;
      width: 13.5rem;
    }
    .region {
      border: 1px solid black;
      border-collapse: collapse;
      display: grid;
      grid-template-columns: 1.5rem 1.5rem 1.5rem;
      grid-template-rows: 1.5rem 1.5rem 1.5rem;
    }
    .cell {
      border: 1px solid gray;
      text-align: center;
    }
    .cell.new {
      color: green;
    }
  </style>
';
$htmlFooter = "</body>\n</html>";
```

Then, to test it out, check out our example instantiation and some code to generate an XML of the puzzle.

```php
// Test it out! This is an "extreme" level difficulty Sudoku puzzle, which is solveable by our dual-algorithmic solver.
$dataStr = "004500003050100008030600700000000100001408096090030000000020057000804200500000000";
// Process this data string into an array of values.
$dataArr = str_split($dataStr);
// Instantiate our puzzle with the data array to represent its initial state.
$puzzle  = new SudokuPuzzle($dataArr);
// If we must solve this puzzle, do so and check for validity.
$isSolved = $puzzle->solve();
// If this puzzle has not been solved, let the user know.
if (!$isSolved) {
  die('<p style="color:red"><strong>ERROR!</strong> Unsolvable '
    . 'puzzle!</p>' . $puzzle->getSolutionHtml());
  }
}

// Begin the SVG output by specifying the SVG content type.
header('Content-Type: image/svg+xml');
// If we must save this SVG as a file for download, use the appropriate headers.
if ($areWeGoingToSaveThisAsAFile) {
  // No cache. Bad cache.
  header('Cache-Control: no-store, no-cache');
  // Generate unique ID based on MD5 hash of data string, ensuring uniqueness.
  $uniqueId = md5($dataStr);
  $puzzle_type = $whichOneDoWeNeed ? 'solution' : 'puzzle';
  // The magic sauce that downloads what we're generating.
  header("Content-Disposition: attachment; "
       . "filename=\"sudoku-puzzle--{$uniqueId}--{$puzzle_type}.svg\"");
}

// Show the puzzle's (or the solution's) SVG code.
echo $puzzle->getPuzzleSvg(/* ... */); // see below for options
```

Here are the three different ways of generating the puzzle's SVG:

```php
// Output the original puzzle, with only the original numbers shown.
echo $puzzle->getPuzzleSvg();
// The default is false, which gets us the base puzzle.

// Output the solution to the puzzle, all numbers populating the grid.
echo $puzzle->getPuzzleSvg(true);

// Output the solution, but with the original numbers a bit lighter,
//  to highlight solved digits in black.
echo $puzzle->getPuzzleSvg(true, true);
```

And that's about it! Updates to come, including additional algorithms as part of the grid passthrough.
