
package sc.fiji.analyzeSkeleton.ita;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Utility methods shared by the test classes.
 *
 * @author Richard Domander (Royal Veterinary College, London)
 * @author Alessandro Felder (Royal Veterinary College, London)
 */
class TestUtil {

	/**
	 * Creates a {@link Graph}, and adds the given edges and vertices to it.
	 * <p>
	 * NB Adds edges as branches of their end points. The connections between the
	 * vertices are defined in the {@link Edge} and {@link Vertex} classes.
	 * </p>
	 *
	 * @see Edge#getV1()
	 * @see Edge#getV2()
	 * @see Vertex#getBranches()
	 * @param edges edges of the graph.
	 * @param vertices vertices of the graph.
	 * @return A graph that contains the vertices and the edge.
	 */
	static Graph createGraph(final Iterable<Edge> edges,
		final Iterable<Vertex> vertices)
	{
		final Graph graph = new Graph();
		edges.forEach(graph::addEdge);
		vertices.forEach(graph::addVertex);
		return graph;
	}
}
