
package sc.fiji.analyzeSkeleton.ita;

import static java.util.Collections.singletonList;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Tests for {@link NJunctions}.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
public class NJunctionsTest {

	@Test
	public void testNJunctionAnglesFrameGraph() {
		final Graph frameGraph = createFrameGraph();

		final List<List<Double>> nJunctions = NJunctions.getNJunctionAngles(
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

		final List<List<Double>> nJunctions = NJunctions.getNJunctionAngles(
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

		final List<List<Double>> nJunctions = NJunctions.getNJunctionAngles(
			vertices);

		assertEquals(2, nJunctions.size());
		for (final List<Double> junctionAngles : nJunctions) {
			assertNotNull(junctionAngles);
			assertTrue(junctionAngles.isEmpty());
		}
	}

	// region -- Helper methods
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
	// endregion
}
