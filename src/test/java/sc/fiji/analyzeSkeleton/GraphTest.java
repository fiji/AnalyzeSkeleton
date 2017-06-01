
package sc.fiji.analyzeSkeleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

/**
 * Tests for {@link Graph}
 *
 * @author Richard Domander
 */
public class GraphTest {

	@Test
	public void testClone() throws Exception {
		final Graph graph = createTestGraph();

		final Graph clone = graph.clone();

		assertTrue(graph != clone);
		final ArrayList<Vertex> vertices = graph.getVertices();
		final ArrayList<Vertex> cloneVertices = clone.getVertices();
		final ArrayList<Edge> edges = graph.getEdges();
		final ArrayList<Edge> cloneEdges = clone.getEdges();
		assertEquals(vertices.size(), cloneVertices.size());
		for (int i = 0; i < vertices.size(); i++) {
			final Vertex vertex = vertices.get(i);
			final Vertex cloneVertex = cloneVertices.get(i);
			assertEquals(vertex.getBranches().size(), cloneVertex.getBranches()
				.size());
			for (int j = 0; j < vertex.getBranches().size(); j++) {
				final Edge branch = vertex.getBranches().get(j);
				final Edge cloneBranch = cloneVertex.getBranches().get(j);
				assertEquals(edges.indexOf(branch), cloneEdges.indexOf(cloneBranch));
			}
		}
		assertTrue("Clone has branches that are not listed in its edges",
			cloneVertices.stream().flatMap(v -> v.getBranches().stream()).allMatch(
				clone.getEdges()::contains));
		assertEquals(cloneEdges.get(0), cloneVertices.get(0).getPredecessor());
		assertEquals(cloneEdges.get(1), cloneVertices.get(1).getPredecessor());
		assertEquals(cloneEdges.get(2), cloneVertices.get(2).getPredecessor());
		assertEquals(edges.size(), cloneEdges.size());
		assertEquals(vertices.indexOf(graph.getRoot()), cloneVertices.indexOf(clone
			.getRoot()));
		for (int i = 0; i < edges.size(); i++) {
			final Edge edge = edges.get(i);
			final Edge cloneEdge = cloneEdges.get(i);
			assertEquals(vertices.indexOf(edge.getV1()), cloneVertices.indexOf(
				cloneEdge.getV1()));
			assertEquals(vertices.indexOf(edge.getV2()), cloneVertices.indexOf(
				cloneEdge.getV2()));
		}
	}

	/**
	 * Creates a graph with 3 vertices and 2 edges, plus a loop edge, and a
	 * parallel edge
	 */
	private Graph createTestGraph() {
		final Graph graph = new Graph();
		final Vertex v1 = new Vertex();
		final Vertex v2 = new Vertex();
		final Vertex v3 = new Vertex();
		final List<Vertex> vertices = Arrays.asList(v1, v2, v3);
		final List<Edge> edges = Arrays.asList(new Edge(v1, v1, null, 0), new Edge(
			v1, v2, null, 0), new Edge(v1, v2, null, 0), new Edge(v2, v3, null, 0));
		v1.setPredecessor(edges.get(0));
		v2.setPredecessor(edges.get(1));
		v3.setPredecessor(edges.get(2));
		edges.forEach(graph::addEdge);
		vertices.forEach(graph::addVertex);
		graph.setRoot(v1);
		return graph;
	}

}
