
package sc.fiji.analyzeSkeleton.ita;

import static java.util.Collections.singletonList;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static sc.fiji.analyzeSkeleton.ita.VertexUtils.groupByValence;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Tests for {@link VertexUtils}.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
public class VertexUtilsTest {

	@Test
	public void testGroupByValenceEmptyCollection() {
		final Map<Integer, List<Vertex>> map = groupByValence(Collections
			.emptyList(), 3, 5);

		assertNotNull(map);
		assertTrue("Empty collection should return an empty map.", map.isEmpty());
	}

	@Test
	public void testGroupByValenceFishGraph() {
		// SETUP
		final int[] expectedCounts = { 2, 2, 1 };
		final int[] expectedValences = { 1, 2, 4 };
		final Graph fishGraph = createDownWardFacingFishGraph();

		// EXECUTE
		final Map<Integer, List<Vertex>> map = groupByValence(fishGraph
			.getVertices(), 0, Integer.MAX_VALUE);

		// VERIFY
		final Set<Integer> keys = map.keySet();
		assertEquals("Wrong number of valences", expectedValences.length, keys
			.size());
		for (int i = 0; i < expectedValences.length; i++) {
			final int expectedValence = expectedValences[i];
			assertTrue("Expected valence missing", map.containsKey(expectedValence));
			final List<Vertex> vertices = map.get(expectedValence);
			assertEquals("Wrong number of vertices with valence " + expectedValence,
				expectedCounts[i], vertices.size());
			assertTrue("A vertex has wrong number of branches (valence)", vertices
				.stream().allMatch(v -> expectedValence == v.getBranches().size()));
		}
	}

	// Tests valence grouping with graph that doesn't have all the valences in the
	// range
	@Test
	public void testGroupByValenceMissingValences() {
		final Graph stickFigureGraph = createStickFigureGraph();

		final Set<Integer> keys = groupByValence(stickFigureGraph.getVertices(), 3,
			5).keySet();

		assertEquals(2, keys.size());
		assertTrue(keys.contains(3));
		assertTrue(keys.contains(5));
	}

	@Test
	public void testGroupByValenceRangeOnFishGraph() {
		final int expectedValence = 2;
		final Graph fishGraph = createDownWardFacingFishGraph();

		final Map<Integer, List<Vertex>> map = groupByValence(fishGraph
			.getVertices(), 2, 3);

		assertEquals("Wrong number of valences", 1, map.keySet().size());
		assertTrue("Expected valence missing", map.containsKey(expectedValence));
	}

	@Test(expected = IllegalArgumentException.class)
	public void testGroupByValenceThrowsIAEIfMinGTMax() {
		groupByValence(singletonList(new Vertex()), 1, 0);
	}

	@Test(expected = IllegalArgumentException.class)
	public void testGroupByValenceThrowsIAEIfMinNegative() {
		groupByValence(singletonList(new Vertex()), -1, 1);
	}

	@Test
	public void testNJunctionAnglesFrameGraph() {
		final Graph frameGraph = createFrameGraph();

		final List<List<Double>> nJunctions = VertexUtils.getNJunctionAngles(
			frameGraph.getVertices());

		assertEquals(4, nJunctions.size());
		final List<Double> junctionAngles = nJunctions.get(0);
		assertEquals(3, junctionAngles.size());
		for (final Double angle : junctionAngles) {
			assertEquals(Math.PI / 2.0, angle, 1e-12);
		}
	}

	@Test
	public void testNJunctionAnglesLonelyVertex() {
		final Vertex lonelyVertex = new Vertex();

		final List<List<Double>> nJunctions = VertexUtils.getNJunctionAngles(
			singletonList(lonelyVertex));

		assertNotNull(nJunctions);
		assertEquals(1, nJunctions.size());
		assertNotNull(nJunctions.get(0));
		assertTrue(nJunctions.get(0).isEmpty());
	}

	@Test
	public void testOneEdge() {
		final List<Vertex> vertices = Arrays.asList(new Vertex(), new Vertex());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(1, 0, 0));
		final Edge edge = new Edge(vertices.get(0), vertices.get(1), null, 1.0);
		// Create a graph, so that edge gets registered as a branch of the vertices
		// etc.
		TestUtil.createGraph(singletonList(edge), vertices);

		final List<List<Double>> nJunctions = VertexUtils.getNJunctionAngles(
			vertices);

		assertEquals(2, nJunctions.size());
		for (final List<Double> junctionAngles : nJunctions) {
			assertNotNull(junctionAngles);
			assertTrue(junctionAngles.isEmpty());
		}
	}

	// region -- Helper methods

	/**
	 * Structure of the downward facing fish graph example
	 *
	 * <pre>
	 * 1   5
	 *  \ /
	 *   2
	 *  / \
	 * 4---3
	 * </pre>
	 */
	private static Graph createDownWardFacingFishGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(5).collect(
			Collectors.toList());
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 0.0), new Edge(vertices.get(1), vertices.get(2), null,
				0.0), new Edge(vertices.get(2), vertices.get(3), null, 0.0), new Edge(
					vertices.get(1), vertices.get(3), null, 0.0), new Edge(vertices.get(
						1), vertices.get(4), null, 0.0));
		return TestUtil.createGraph(edges, vertices);
	}

	/**
	 * Structure of the frame graph example ("frame" as in "coordinate frame")
	 *
	 * <pre>
	 * z
	 * |  y
	 * | /
	 * |/
	 * o-----x
	 * </pre>
	 */
	private static Graph createFrameGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(4).collect(
			Collectors.toList());

		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(1, 0, 0));
		vertices.get(2).addPoint(new Point(0, 1, 0));
		vertices.get(3).addPoint(new Point(0, 0, 1));

		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 1.0), new Edge(vertices.get(0), vertices.get(2), null,
				1.0), new Edge(vertices.get(0), vertices.get(3), null, 1.0));

		return TestUtil.createGraph(edges, vertices);
	}

	/**
	 * Structure of the stick figure graph example
	 *
	 * <pre>
	 * 0     1
	 *  \   /
	 *   \ /
	 * 6--2--7
	 *    |
	 *    |
	 *    |
	 *    |
	 *    3
	 *   / \
	 *  /   \
	 * 4     5
	 * </pre>
	 */
	private static Graph createStickFigureGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(8).collect(
			Collectors.toList());

		vertices.get(0).addPoint(new Point(-2, 4, 0));
		vertices.get(1).addPoint(new Point(2, 4, 0));
		vertices.get(2).addPoint(new Point(0, 3, 0));
		vertices.get(3).addPoint(new Point(0, -3, 0));
		vertices.get(4).addPoint(new Point(-2, -4, 0));
		vertices.get(5).addPoint(new Point(-2, -4, 0));
		vertices.get(5).addPoint(new Point(2, 3, 0));
		vertices.get(5).addPoint(new Point(2, 3, 0));

		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(2), null, Math.sqrt(5.0)), new Edge(vertices.get(1), vertices.get(2),
				null, Math.sqrt(5.0)), new Edge(vertices.get(2), vertices.get(3), null,
					6.0), new Edge(vertices.get(3), vertices.get(4), null, Math.sqrt(
						5.0)), new Edge(vertices.get(3), vertices.get(5), null, Math.sqrt(
							5.0)), new Edge(vertices.get(6), vertices.get(2), null, 2.0),
			new Edge(vertices.get(7), vertices.get(2), null, 2.0));

		return TestUtil.createGraph(edges, vertices);
	}
	// endregion
}
