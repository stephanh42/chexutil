import unittest
import chexutil
import pickle
import math

class HexMap(object):
    def __init__(self, str):
        self.source = str
        player = None
        target = None
        tiles = {}
        line_lengths = []
        lights = []
        for y, line in enumerate(str.split("\n")):
            line_lengths.append(len(line))
            for x, ch in enumerate(line):
                if ch.isspace():
                    continue
                position = chexutil.Hex(x, y)
                if ch == '@':
                    player = position
                    ch = '.'
                elif ch == '%':
                    ch = '.'
                elif ch == 'T':
                    target = position
                elif ch == '*':
                    lights.append(position)
                tiles[position] = ch
        self.player = player
        self.target = target
        self.tiles = tiles
        self.line_lengths = line_lengths
        self.lights = lights

    def is_transparent(self, pos):
        return self.tiles.get(pos, '#') != '#'

    def is_passable(self, pos):
        return self.tiles.get(pos, '#') not in "#~"

    def field_of_view(self, max_distance):
        return self.player.field_of_view(transparent=self.is_transparent, max_distance=max_distance)

    def get_map(self, max_distance=None, path=frozenset()):
        if max_distance is not None:
            fov = self.field_of_view(max_distance)
            if self.lights:
                light_fov = {}
                for light in self.lights:
                    light.field_of_view(transparent=self.is_transparent, max_distance=max_distance, visible=light_fov)
            else:
                light_fov = fov
        else:
            fov = light_fov = None
        lines = []
        for y, line_length in enumerate(self.line_lengths):
            line = []
            for x in range(line_length):
                if (x + y) % 2 != 0:
                    ch = ' '
                else:
                    pos = chexutil.Hex(x, y)
                    if pos == self.player:
                        ch = '@'
                    elif fov is None or (fov.get(pos, 0) & light_fov.get(pos, 0)):
                        if pos in path:
                            ch = '%'
                        else:
                            ch = self.tiles.get(pos, ' ')
                    else:
                        ch = ' '
                line.append(ch)
            lines.append("".join(line).rstrip())
        return "\n".join(lines)


testmap1 = HexMap("""
 # # # # # # # # #
# # # # # # # # # #
 # . # # # # # # #
# # . # # # # # # #
 # # # . . . . # #
# # # . @ # # . # #
 # # ~ % . # # . #
# # ~ # % # # # # #
 # # T % # # # # #
# # # # # # # # # #
""")

testmap1_out = """



      # # # #
     # . . .
    # . @ #
   # ~ . . #
  # ~ # . #
   #   . #
      # #
"""

testmap2 = HexMap("""
 # # # # # # # # #
# # # # # # # # # #
 # . # # # * # # #
# # . # # . # # # #
 # # # . % % % # #
# # . . @ # # % # #
 # . # . # # # T #
# # # # . . # # # #
 # # . . . . # # #
# # # # # # # # # #
""")

testmap2_out = """

          # #
           * #
          .
       . . .
      . @
     # . #
      # .


"""


class TestHex(unittest.TestCase):

    def test_is_valid(self):
        chexutil.Hex(-1, -3)
        self.assertRaises(chexutil.InvalidHex, chexutil.Hex, 2, -3)

    def test_rotations(self):
        hex = chexutil.Hex(2, 0)
        self.assertEqual([r(hex) for r in chexutil.Hex.rotations], chexutil.origin.neighbours())

    def test_add(self):
        self.assertEqual(chexutil.Hex(2, 4) + chexutil.Hex(4, 6), chexutil.Hex(6, 10))

    def test_sub(self):
        self.assertEqual(chexutil.Hex(2, 4) - chexutil.Hex(3, 7), chexutil.Hex(-1, -3))

    def test_neg(self):
        self.assertEqual(-chexutil.Hex(2, 4), chexutil.Hex(-2, -4))
        self.assertEqual(-chexutil.origin, chexutil.origin)

    def test_neighbours(self):
        origin = chexutil.origin
        nb = origin.neighbours()
        self.assertEqual(len(nb), 6)
        for h in nb:
            self.assertTrue((-h) in nb)
            self.assertTrue(origin in h.neighbours())

    def test_distance(self):
        h = chexutil.Hex(-1, -3)
        for nb in h.neighbours():
            self.assertEqual(nb.distance(nb), 0)
            self.assertEqual(nb.distance(h), 1)
            self.assertEqual(max(nb2.distance(h) for nb2 in nb.neighbours()), 2)

    def test_rotate_left(self):
        origin = chexutil.origin
        neighbours = origin.neighbours()
        nb = neighbours[0]
        neighbours2 = []
        for i in range(6):
            neighbours2.append(nb)
            nb = nb.rotate_left()
        self.assertEqual(neighbours, neighbours2)

    def test_rotate_right(self):
        origin = chexutil.origin
        for nb in chexutil.origin.neighbours():
            self.assertEqual(nb.rotate_left().rotate_right(), nb)

class TestHexGrid(unittest.TestCase):
    def test_height(self):
        self.assertEqual(chexutil.HexGrid(32).height, 18)

    def test_corners(self):
        self.assertEqual(chexutil.HexGrid(32).corners(chexutil.Hex(1, 1)), 
                 [(64, 72), (32, 90), (0, 72), (0, 36), (32, 18), (64, 36)])

    def test_hex_at_coordinates(self):
        hg = chexutil.HexGrid(32)
        data = [
                ((0, 0), chexutil.Hex(0, 0)),
                ((33, 16), chexutil.Hex(2, 0)),
                ((30, 20), chexutil.Hex(1, 1))
                ]
        for fx in (-1, 1):
            for fy in (-1, 1):
                for pixel, hex in data:
                    x, y = pixel
                    pixel = (fx*x, fy*y)
                    x, y = hex
                    hex = chexutil.Hex(fx*x, fy*y)
                    self.assertEqual(hg.hex_at_coordinate(*pixel), hex, pixel)

    def test_center(self):
        hg = chexutil.HexGrid(32)
        self.assertEqual(hg.center(chexutil.Hex(1, 1)), (32, 54))

    def test_bounding_box(self):
        hg = chexutil.HexGrid(32)
        self.assertEqual(hg.bounding_box(chexutil.Hex(0,2)),
                chexutil.Rectangle(-32, 72, 64, 72))

    def test_hexes_in_rectangle(self):
        hg = chexutil.HexGrid(32)
        self.assertEqual(
                list(hg.hexes_in_rectangle(hg.bounding_box(chexutil.origin))),
                [chexutil.Hex(-1, -1), chexutil.Hex(1, -1),
                 chexutil.Hex(-2, 0), chexutil.Hex(0, 0),
                 chexutil.Hex(-1, 1), chexutil.Hex(1, 1)]
                )

    def test_hexrange(self):
        for x1, x2, y1, y2 in [
                (1, 5, 2, 4),
                (1, 3, 2, 5),
                (2, 3, 1, 5),
                (1, 5, 2, 3),
                (2, 3, 0, 5),
                (0, 5, 2, 3),
                (0, 3, 2, 2)]:
            r = chexutil.HexRange(x1, x2, y1, y2)
            l = [chexutil.Hex(x, y) for y in range(y1, y2) for x in range(x1, x2) if (x + y) % 2 == 0]
            self.assertEqual(list(r), l)

    def test_hex_factor(self):
        self.assertEqual(chexutil.hex_factor, math.sqrt(1.0/3.0))

class TestFov(unittest.TestCase):
    def test_fov1(self):
        self.assertEqual(testmap1.get_map(10), testmap1_out)

    def test_fov2(self):
        self.assertEqual(testmap2.get_map(10), testmap2_out)

class TestPathFinding(unittest.TestCase):
    def test_path1(self):
        path = frozenset(testmap1.player.find_path(testmap1.target, testmap1.is_passable)[:-1])
        self.assertEqual(testmap1.get_map(path=path), testmap1.source)

    def test_path2(self):
        path = frozenset(testmap2.player.find_path(testmap2.target, testmap2.is_passable)[:-1])
        self.assertEqual(testmap2.get_map(path=path), testmap2.source)

class TextRect(unittest.TestCase):
    def test_rect1(self):
        self.assertEqual(chexutil.Rectangle(4, 5).translated(1, 2), chexutil.Rectangle(1, 2, 4, 5))

    def test_union(self):
        R = chexutil.Rectangle
        self.assertEqual(R(1, 2, 3, 4) | R(4, 5, 1, 2), R(1, 2, 4, 5))

    def test_intersect(self):
        R = chexutil.Rectangle
        self.assertEqual(R(1, 2, 3, 4) & R(4, 5, 1, 2), R(0, 0, 0, 0))
        self.assertEqual(R(1, 2, 5, 6) & R(4, 5, 1, 2), R(4, 5, 1, 2))

class TestPickle(unittest.TestCase):
    def test_pickle(self):
        for obj in [chexutil.Hex(1, 3), chexutil.Rectangle(1, 2, 3, 4)]:
            self.assertEqual(obj, pickle.loads(pickle.dumps(obj)))

if __name__ == '__main__':
    unittest.main()
