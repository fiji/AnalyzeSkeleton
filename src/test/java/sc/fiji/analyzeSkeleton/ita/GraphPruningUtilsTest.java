
package sc.fiji.analyzeSkeleton.ita;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.findClusters;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.findEdgesWithOneEndInCluster;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.getClusterCentre;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.isShort;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.pruneShortEdges;
import static sc.fiji.analyzeSkeleton.ita.GraphPruningUtils.removeParallelEdges;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.junit.Test;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Tests for {@link GraphPruningUtils}.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
public class GraphPruningUtilsTest {

	private static final double[] ISOTROPIC = { 1.0, 1.0, 1.0 };

	@Test
	public void testFindClustersCentresHaveOnePoint() {
		// SETUP
		final Graph graph = createDumbbellGraph();
		final List<Set<Vertex>> clusters = findClusters(graph, 2.01);

		// EXECUTE
		final List<Vertex> centres = clusters.stream().map(
			GraphPruningUtils::getClusterCentre).collect(Collectors.toList());

		// VERIFY
		assertTrue(centres.stream().allMatch(c -> c.getPoints().size() == 1));
	}

	// Should result in a graph with two clusters of 3 connected by a single edge
	@Test
	public void testFindClustersDumbbell() {
		// SETUP
		final Graph dumbbellGraph = createDumbbellGraph();

		// EXECUTE
		final List<Set<Vertex>> clusters = findClusters(dumbbellGraph, 2.01);

		// VERIFY CLUSTERS
		assertNotNull(clusters);
		assertEquals(2, clusters.size());
		final Set<Vertex> firstCluster = clusters.get(0);
		assertEquals(3, firstCluster.size());
		final Set<Vertex> secondCluster = clusters.get(1);
		assertEquals(3, secondCluster.size());
		// VERIFY CLUSTER VERTICES
		final ArrayList<Vertex> vertices = dumbbellGraph.getVertices();
		assertTrue("Clusters share vertices", vertices.stream().filter(
			v -> !firstCluster.contains(v)).allMatch(secondCluster::contains));
		assertTrue("Clusters share vertices", vertices.stream().filter(
			v -> !secondCluster.contains(v)).allMatch(firstCluster::contains));
		// VERIFY CLUSTER EDGES
		final List<Edge> firstClusterEdges = firstCluster.stream().flatMap(v -> v
			.getBranches().stream()).filter(e -> isShort(e, 2.01)).collect(Collectors
				.toList());
		final List<Edge> secondClusterEdges = secondCluster.stream().flatMap(v -> v
			.getBranches().stream()).filter(e -> isShort(e, 2.01)).collect(Collectors
				.toList());
		final List<Edge> allShortEdges = dumbbellGraph.getEdges().subList(0, 6);
		assertTrue("There are short edges that are not contained in either cluster",
			allShortEdges.stream().filter(e -> !firstClusterEdges.contains(e))
				.allMatch(secondClusterEdges::contains));
		final Edge longEdge = dumbbellGraph.getEdges().get(6);
		assertEquals("There should be one edge connecting the clusters", 1,
			firstCluster.stream().flatMap(v -> v.getBranches().stream()).filter(
				longEdge::equals).count());
		assertTrue("There are short edges that are not contained in either cluster",
			allShortEdges.stream().filter(e -> !secondClusterEdges.contains(e))
				.allMatch(firstClusterEdges::contains));
		assertEquals("There should be one edge connecting the clusters", 1,
			secondCluster.stream().flatMap(v -> v.getBranches().stream()).filter(
				longEdge::equals).count());
	}

	@Test
	public void testFindClustersKite() {
		// SETUP
		final Graph kiteGraph = createKiteGraph();
		final ArrayList<Vertex> vertices = kiteGraph.getVertices();
		final ArrayList<Edge> edges = kiteGraph.getEdges();

		// EXECUTE
		final List<Set<Vertex>> clusters = findClusters(kiteGraph, 1.01);

		// VERIFY
		assertNotNull(clusters);
		assertEquals(1, clusters.size());
		final Set<Vertex> cluster = clusters.get(0);
		assertEquals(3, cluster.size());
		assertTrue("Clustering created unexpected vertices", vertices.stream()
			.limit(3).allMatch(cluster::contains));
		final Set<Edge> clusterEdgeSet = cluster.stream().flatMap(v -> v
			.getBranches().stream()).collect(Collectors.toSet());
		assertEquals(4, clusterEdgeSet.size());
		assertTrue("Clustering created unexpected edges", edges.stream().limit(4)
			.allMatch(clusterEdgeSet::contains));
	}

	@Test
	public void testFindClustersSquareCluster() {
		// SETUP
		final Graph squareWithDiagAndLooseEnds = createTriangleWithSquareCluster();

		// EXECUTE
		final List<Set<Vertex>> clusters = findClusters(squareWithDiagAndLooseEnds,
			2.01);

		// VERIFY
		assertEquals(1, clusters.size());
		final Set<Vertex> cluster = clusters.get(0);
		assertEquals(4, cluster.size());
		assertEquals(7, cluster.stream().map(Vertex::getBranches).flatMap(
			List::stream).distinct().count());
	}

	@Test
	public void testFindEdgesWithOneEndInClusterSquareCluster() {
		// SETUP
		final Graph squareWithDiagAndLooseEnds = createTriangleWithSquareCluster();
		final Set<Vertex> cluster = findClusters(squareWithDiagAndLooseEnds, 2.01)
			.get(0);

		// EXECUTE
		final Set<Edge> outerEdges = findEdgesWithOneEndInCluster(cluster);

		// VERIFY
		assertEquals(
			"Cluster has wrong number of edges connecting it to the rest of the graph",
			2, outerEdges.size());
	}

	@Test
	public void testGetClusterCentreCentroidCoordinates() {
		// SETUP
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(3).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(2, 0, 0));
		vertices.get(1).addPoint(new Point(4, 0, 0));
		vertices.get(2).addPoint(new Point(3, 6, 0));
		final Set<Vertex> cluster = new HashSet<>(vertices);

		// EXECUTE
		final Vertex clusterCentre = getClusterCentre(cluster);

		// VERIFY
		assertTrue("Centroid has wrong coordinates", clusterCentre.isVertexPoint(
			new Point(3, 2, 0)));
	}

	@Test
	public void testGetClusterCentreEmptyCluster() {
		// EXECUTE
		final ArrayList<Point> points = getClusterCentre(Collections.emptySet())
			.getPoints();

		// VERIFY
		assertEquals(1, points.size());
		assertEquals(Integer.MAX_VALUE, points.get(0).x);
		assertEquals(Integer.MAX_VALUE, points.get(0).y);
		assertEquals(Integer.MAX_VALUE, points.get(0).z);
	}

	@Test
	public void testGetClusterCentreIsSinglePoint() {
		// SETUP
		final Graph graph = createDumbbellGraph();
		final List<Set<Vertex>> clusters = findClusters(graph, 2.01);

		// EXECUTE
		final List<Vertex> centres = clusters.stream().map(
			GraphPruningUtils::getClusterCentre).collect(Collectors.toList());

		// VERIFY
		assertTrue(centres.stream().allMatch(c -> c.getPoints().size() == 1));
	}

	@Test
	public void testPruneShortEdgesAnisotropicVoxels() {
		// SETUP
		final Graph sailGraph = createSailGraph();
		final double[] calibration = { 2.0, 5.0, 3.0 };

		// EXECUTE
		final Graph cleanSailGraph = pruneShortEdges(sailGraph, 0, false, true,
			calibration);

		// VERIFY
		assertEquals(4, cleanSailGraph.getEdges().size());
		assertTrue(cleanSailGraph.getEdges().stream().anyMatch(e -> e
			.getLength() == 4.0));
		assertTrue(cleanSailGraph.getEdges().stream().anyMatch(e -> e
			.getLength() == 15.0));
		assertTrue(cleanSailGraph.getEdges().stream().anyMatch(e -> e
			.getLength() == Math.sqrt(241)));
		assertTrue(cleanSailGraph.getEdges().stream().anyMatch(e -> e
			.getLength() == 5.0));
	}

	// Tests that pruning works differently, when it operates on vertex pairs
	// instead of clusters of vertices
	@Test
	public void testPruneShortEdgesClusteredFalse() {
		// SETUP
		final Graph segmentGraph = createLineGraph();

		// EXECUTE
		final Graph cleanSegmentGraph = GraphPruningUtils.pruneShortEdges(
			segmentGraph, 4.01, false, false, ISOTROPIC);

		// VERIFY
		assertEquals(3, cleanSegmentGraph.getVertices().size());
		assertEquals(2, cleanSegmentGraph.getEdges().size());
		assertEquals(1, cleanSegmentGraph.getEdges().stream().filter(e -> e
			.getLength() == 17.0).count());
		assertEquals(1, cleanSegmentGraph.getEdges().stream().filter(e -> e
			.getLength() == 15.0).count());
	}

	// Tests that edges created by the pruning have linear length, i.e. the
	// euclidean distance of the centroids of their end points.
	@Test
	public void testPruneShortEdgesEdgeLength() {
		// SETUP
		final Graph arch = createSlabArchGraph();

		// EXECUTE
		final Graph cleanArch = pruneShortEdges(arch, 0, false, true);

		// VERIFY
		assertEquals(5, cleanArch.getEdges().get(0).getLength(), 1e-12);
	}

	@Test
	public void testPruneShortEdgesEmptyGraph() {
		// SETUP
		final Graph emptyGraph = new Graph();

		// EXECUTE
		final Graph cleanEmptyGraph = pruneShortEdges(emptyGraph, 2.01, false, true,
			new double[] {});

		// VERIFY
		assertNotNull(cleanEmptyGraph);
		assertTrue(cleanEmptyGraph.getVertices().isEmpty());
		assertTrue(cleanEmptyGraph.getEdges().isEmpty());
	}

	@Test
	public void testPruneShortEdgesInputNotMutated() {
		// SETUP
		final Graph graph = createSailGraph();
		final Graph cloneGraph = graph.clone();

		// EXECUTE
		pruneShortEdges(graph, Double.POSITIVE_INFINITY, false, true, ISOTROPIC);

		// VERIFY
		assertGraphEquals(graph, cloneGraph);
	}

	// Tests that pruning a graph more than once creates a different result
	@Test
	public void testPruneShortEdgesIterativePruning() {
		// SETUP
		final Graph doorknob = createDoorknobGraph();

		// EXECUTE
		final Graph cleanedOnce = pruneShortEdges(doorknob, 2.01, false, true,
			ISOTROPIC);
		final Graph cleanedTwice = pruneShortEdges(doorknob, 2.01, true, true,
			ISOTROPIC);

		// VERIFY
		assertEquals(4, cleanedOnce.getVertices().size());
		assertEquals(3, cleanedOnce.getEdges().size());
		assertEquals(3, cleanedTwice.getVertices().size());
		assertEquals(2, cleanedTwice.getEdges().size());
	}

	@Test
	public void testPruneShortEdgesLonelyVertexReturnsLonelyVertex() {
		// SETUP
		final Vertex vertex = new Vertex();
		final Point randomNonZero = new Point(7, 1, 3);
		vertex.addPoint(randomNonZero);
		final Graph oneVertexGraph = new Graph();
		oneVertexGraph.addVertex(vertex);

		// EXECUTE
		final Graph cleanOneVertexGraph = pruneShortEdges(oneVertexGraph,
			Double.POSITIVE_INFINITY, false, true, new double[] {});

		// VERIFY
		assertEquals(1, cleanOneVertexGraph.getVertices().size());
		assertEquals(1, cleanOneVertexGraph.getVertices().get(0).getPoints()
			.size());
		assertEquals(randomNonZero, cleanOneVertexGraph.getVertices().get(0)
			.getPoints().get(0));
	}

	@Test
	public void testPruneShortEdgesNoParallelEdgesInOutput() {
		// SETUP
		final Graph graph = createTriangleWithSquareClusterAndArtefacts();

		// EXECUTE
		final Graph result = pruneShortEdges(graph, Double.POSITIVE_INFINITY, false,
			true, ISOTROPIC);

		// VERIFY
		final int size = result.getEdges().size();
		removeParallelEdges(result);
		final int secondPassSize = result.getEdges().size();
		assertEquals(size, secondPassSize);
	}

	@Test
	public void testPruneShortEdgesRemovesLoops() {
		// SETUP
		final Graph loopGraph = createLoopGraph();

		// EXECUTE
		final Graph cleanLoopGraph = GraphPruningUtils.pruneShortEdges(loopGraph,
			0.0, false, true, ISOTROPIC);

		// VERIFY
		assertEquals(3, cleanLoopGraph.getEdges().size());
		assertTrue(cleanLoopGraph.getEdges().stream().filter(e -> e.getV1() == null)
			.noneMatch(e -> e.getV1() == e.getV2()));
	}

	@Test
	public void testPruneShortEdgesSail() {
		// SETUP
		final Graph sailGraph = createSailGraph();

		// EXECUTE
		final Graph cleanSailGraph = pruneShortEdges(sailGraph, 1.01, false, true,
			ISOTROPIC);

		// VERIFY
		assertEquals(3, cleanSailGraph.getEdges().size());
		assertEquals(3, cleanSailGraph.getVertices().size());
		assertEquals(1, cleanSailGraph.getVertices().stream().flatMap(v -> v
			.getPoints().stream()).filter(p -> p.x == 0 && p.y == 0 && p.z == 0)
			.count());
		assertEquals(1, cleanSailGraph.getVertices().stream().flatMap(v -> v
			.getPoints().stream()).filter(p -> p.x == 0 && p.y == 3 && p.z == 0)
			.count());
		assertEquals(1, cleanSailGraph.getVertices().stream().flatMap(v -> v
			.getPoints().stream()).filter(p -> p.x == 2 && p.y == 0 && p.z == 0)
			.count());
	}

	@Test
	public void testPruneShortEdgesTriangleWithSquareCluster() {
		// SETUP
		final Graph squareWithDiagAndLooseEnds = createTriangleWithSquareCluster();
		final List<Point> expectedPoints = Arrays.asList(new Point(0, 0, 0),
			new Point(5, 4, 0), new Point(-4, -5, 0), new Point(5, -5, 0));

		// EXECUTE
		final Graph cleaned = pruneShortEdges(squareWithDiagAndLooseEnds, 2.01,
			false, true, ISOTROPIC);

		// VERIFY
		assertNotNull(cleaned);
		assertEquals(4, cleaned.getVertices().size());
		assertEquals(4, cleaned.getEdges().size());
		expectedPoints.forEach(p -> {
			assertEquals("There should be exactly one vertex with the given point", 1,
				cleaned.getVertices().stream().filter(v -> hasPointOnce(v, p)).count());
			assertEquals("There should be exactly 2 edges with the given point", 2,
				cleaned.getEdges().stream().filter(e -> hasPoint(e, p)).count());
		});
		final double expectedLength1 = Math.sqrt(41.0);
		final double expectedLength2 = 9.0;
		assertEquals(2, cleaned.getEdges().stream().filter(e -> e
			.getLength() == expectedLength1).count());
		assertEquals(2, cleaned.getEdges().stream().filter(e -> e
			.getLength() == expectedLength2).count());
	}

	@Test
	public void testRemoveParallelEdges() {
		// SETUP
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(2).collect(
			Collectors.toList());
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 0), new Edge(vertices.get(1), vertices.get(0), null, 0),
			new Edge(vertices.get(0), vertices.get(1), null, 0));
		final Graph graph = createGraph(edges, vertices);

		// EXECUTE
		removeParallelEdges(graph);

		// VERIFY
		assertEquals(1, graph.getEdges().size());
		assertTrue(graph.getVertices().stream().map(v -> v.getBranches().size())
			.allMatch(b -> b == 1));
	}

	private static void assertEdgeEquals(final Edge expected, final Edge actual,
		final List<Vertex> expectedVertices, final List<Vertex> actualVertices)
	{
		assertEquals(expected.getColor(), actual.getColor(), 1e-12);
		assertEquals(expected.getColor3rd(), actual.getColor3rd(), 1e-12);
		assertEquals(expected.getLength(), actual.getLength(), 1e-12);
		assertEquals(expected.getLength_ra(), actual.getLength_ra(), 1e-12);
		assertEquals(expected.getType(), actual.getType(), 1e-12);
		assertEquals(expectedVertices.indexOf(expected.getV1()), actualVertices
			.indexOf(actual.getV1()));
		assertEquals(expectedVertices.indexOf(expected.getV2()), actualVertices
			.indexOf(actual.getV2()));
	}

	// region -- Helper regions --

	private static void assertGraphEquals(final Graph expected,
		final Graph actual)
	{
		final ArrayList<Vertex> vertices = expected.getVertices();
		final ArrayList<Vertex> cloneVertices = actual.getVertices();
		assertEquals(vertices.size(), cloneVertices.size());
		IntStream.range(0, vertices.size()).forEach(i -> {
			final Vertex vertex = vertices.get(i);
			final Vertex clonedVertex = cloneVertices.get(i);
			assertVertexEquals(vertex, clonedVertex, expected.getEdges(), actual
				.getEdges());
		});
		// VERIFY EDGES
		final ArrayList<Edge> clonedEdges = actual.getEdges();
		assertTrue("Vertices have branches that are not listed in the Edge list",
			cloneVertices.stream().flatMap(v -> v.getBranches().stream()).allMatch(
				clonedEdges::contains));
		final ArrayList<Edge> edges = expected.getEdges();
		assertEquals("Graph has wrong number of edges", edges.size(), clonedEdges
			.size());
		IntStream.range(0, edges.size()).forEach(i -> {
			final Edge edge = edges.get(i);
			final Edge clonedEdge = clonedEdges.get(i);
			assertEdgeEquals(edge, clonedEdge, expected.getVertices(), actual
				.getVertices());
		});
	}

	private static void assertPointsEquals(final List<Point> expected,
		final List<Point> actual)
	{
		assertEquals("Cloned vertex has wrong number of points", expected.size(),
			actual.size());
		IntStream.range(0, expected.size()).forEach(j -> {
			final Point point = expected.get(j);
			final Point clonedPoint = actual.get(j);
			assertEquals("Cloned vertex has a point with bad coordinates", point
				.toString(), clonedPoint.toString());
		});
	}

	private static void assertVertexEquals(final Vertex expected,
		final Vertex actual, final List<Edge> expectedEdges,
		final List<Edge> actualEdges)
	{
		assertEquals(expected.getVisitOrder(), actual.getVisitOrder());
		assertEquals(expected.isVisited(), actual.isVisited());
		assertEquals(expectedEdges.indexOf(expected.getPredecessor()), actualEdges
			.indexOf(actual.getPredecessor()));
		assertEquals("Vertex has wrong number of branches", expected.getBranches()
			.size(), actual.getBranches().size());
		assertPointsEquals(expected.getPoints(), actual.getPoints());
	}

	/**
	 * Creates a {@link Graph} shaped like a door knob.
	 *
	 * <pre>
	 *   4
	 *   |
	 *   |
	 *   |
	 * 1 |
	 * |\2--6
	 * |
	 * |/3--7
	 * 0 |
	 *   |
	 *   |
	 *   |
	 *   5
	 * </pre>
	 */
	private static Graph createDoorknobGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(8).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(0, 3, 0));
		vertices.get(2).addPoint(new Point(1, 2, 0));
		vertices.get(3).addPoint(new Point(1, 1, 0));
		vertices.get(4).addPoint(new Point(1, 6, 0));
		vertices.get(5).addPoint(new Point(1, -3, 0));
		vertices.get(6).addPoint(new Point(3, 2, 0));
		vertices.get(7).addPoint(new Point(3, 1, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 3.0), new Edge(vertices.get(0), vertices.get(3), null, Math
				.sqrt(2.0)), new Edge(vertices.get(1), vertices.get(2), null, Math.sqrt(
					2.0)), new Edge(vertices.get(2), vertices.get(4), null, Math.sqrt(
						4.0)), new Edge(vertices.get(3), vertices.get(5), null, Math.sqrt(
							4.0)), new Edge(vertices.get(2), vertices.get(6), null, 2.0),
			new Edge(vertices.get(3), vertices.get(7), null, 2.0));
		return createGraph(edges, vertices);
	}

	/**
	 * Creates a {@link Graph} shaped like a dumbbell.
	 *
	 * <pre>
	 * 2      5
	 *  \    /
	 *   1--3
	 *  /    \
	 * 0      4
	 * </pre>
	 */
	private static Graph createDumbbellGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(6).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, -1, 0));
		vertices.get(1).addPoint(new Point(1, 0, 0));
		vertices.get(2).addPoint(new Point(0, 1, 0));
		vertices.get(3).addPoint(new Point(4, 0, 0));
		vertices.get(4).addPoint(new Point(5, -1, 0));
		vertices.get(5).addPoint(new Point(5, 1, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, Math.sqrt(2.0)), new Edge(vertices.get(0), vertices.get(2),
				null, 2.0), new Edge(vertices.get(1), vertices.get(2), null, Math.sqrt(
					2.0)), new Edge(vertices.get(3), vertices.get(4), null, Math.sqrt(
						2.0)), new Edge(vertices.get(4), vertices.get(5), null, 2.0),
			new Edge(vertices.get(5), vertices.get(3), null, Math.sqrt(2.0)),
			new Edge(vertices.get(1), vertices.get(3), null, 3.0));

		return createGraph(edges, vertices);
	}

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
	private static Graph createGraph(final Iterable<Edge> edges,
		final Iterable<Vertex> vertices)
	{
		final Graph graph = new Graph();
		edges.forEach(graph::addEdge);
		vertices.forEach(graph::addVertex);
		return graph;
	}

	/**
	 * Creates a {@link Graph} in the shape of a kite.
	 *
	 * <pre>
	 *       ______4
	 *     _/     _/
	 *    /      /
	 *  _/     _/
	 * 2     _3
	 * |   _/
	 * 0--1
	 * </pre>
	 */
	private static Graph createKiteGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(5).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(1, 0, 0));
		vertices.get(2).addPoint(new Point(0, 1, 0));
		vertices.get(3).addPoint(new Point(2, 2, 0));
		vertices.get(4).addPoint(new Point(5, 5, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 1.0), new Edge(vertices.get(0), vertices.get(2), null,
				1.0), new Edge(vertices.get(1), vertices.get(3), null, Math.sqrt(3.0)),
			new Edge(vertices.get(2), vertices.get(3), null, Math.sqrt(3.0)),
			new Edge(vertices.get(3), vertices.get(4), null, 3 * Math.sqrt(2.0)));
		return createGraph(edges, vertices);
	}

	/**
	 * Creates a {@link Graph} with 5 vertices in a straight line.
	 *
	 * <pre>
	 * 0---1-2-3---4
	 * </pre>
	 */
	private static Graph createLineGraph() {

		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(5).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(-16, 0, 0));
		vertices.get(1).addPoint(new Point(-4, 0, 0));
		vertices.get(2).addPoint(new Point(0, 0, 0));
		vertices.get(3).addPoint(new Point(4, 0, 0));
		vertices.get(4).addPoint(new Point(16, 0, 0));

		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 12.0), new Edge(vertices.get(1), vertices.get(2), null,
				4.0), new Edge(vertices.get(2), vertices.get(3), null, 4.0), new Edge(
					vertices.get(3), vertices.get(4), null, 12.0));

		return createGraph(edges, vertices);
	}

	/**
	 * Creates a {@link Graph} with a loop in it.
	 * <p>
	 * "o" denotes a zero length loop edge
	 * </p>
	 *
	 * <pre>
	 * 	 o
	 *   0
	 *  / \
	 * 1---2
	 * </pre>
	 */
	private static Graph createLoopGraph() {

		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(3).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(-1, -1, 0));
		vertices.get(2).addPoint(new Point(1, -1, 0));

		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(0), null, 0.0), new Edge(vertices.get(0), vertices.get(1), null,
				1.0), new Edge(vertices.get(0), vertices.get(2), null, 1.0), new Edge(
					vertices.get(1), vertices.get(2), null, 2.0));

		return createGraph(edges, vertices);
	}

	/**
	 * Creates a {@link Graph} shaped like a sail.
	 * <p>
	 * The graph can be used to test pruning of "dead ends". Edge (0)-(3) is
	 * considered a dead end if tolerance is set &ge; 1.
	 * </p>
	 *
	 * <pre>
	 * 2
	 * |\
	 * | |
	 * | \
	 * 0--1
	 * |
	 * 3
	 * </pre>
	 */
	private static Graph createSailGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(4).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(2, 0, 0));
		vertices.get(2).addPoint(new Point(0, 3, 0));
		vertices.get(3).addPoint(new Point(0, -1, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 2.0), new Edge(vertices.get(0), vertices.get(2), null,
				3.0), new Edge(vertices.get(1), vertices.get(2), null, Math.sqrt(13.0)),
			new Edge(vertices.get(0), vertices.get(3), null, 1.0));
		return createGraph(edges, vertices);
	}

	/**
	 * Creates a {@link Graph} with an arched edge.
	 * <p>
	 * v: vertex point s: slab point of single edge in graph
	 * </p>
	 *
	 * <pre>
	 *   ss
	 * 	s  s
	 * v	v
	 * </pre>
	 */
	private static Graph createSlabArchGraph() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(2).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(0, 0, 0));
		vertices.get(1).addPoint(new Point(5, 0, 0));
		final List<Point> slabPoints = Arrays.asList(new Point(1, 1, 0), new Point(
			2, 2, 0), new Point(3, 2, 0), new Point(4, 1, 0));
		final ArrayList<Point> slabs = new ArrayList<>(slabPoints);
		final Edge edge = new Edge(vertices.get(0), vertices.get(1), slabs, 4 * Math
			.sqrt(2.0) + 1);
		return createGraph(Collections.singleton(edge), vertices);
	}

	/**
	 * Creates a triangle {@link Graph} with square cluster
	 *
	 * <pre>
	 *             4
	 *           _/|
	 *     3--2_/  |
	 *     |\_|    |
	 *    _0--1    |
	 *  _/         |
	 * /           |
	 *5------------6
	 * </pre>
	 */
	private static Graph createTriangleWithSquareCluster() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(7).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(-1, -1, 0));
		vertices.get(1).addPoint(new Point(-1, 1, 0));
		vertices.get(2).addPoint(new Point(1, 1, 0));
		vertices.get(3).addPoint(new Point(1, -1, 0));
		vertices.get(4).addPoint(new Point(5, 4, 0));
		vertices.get(5).addPoint(new Point(-4, -5, 0));
		vertices.get(6).addPoint(new Point(5, -5, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 2.0), new Edge(vertices.get(1), vertices.get(2), null,
				2.0), new Edge(vertices.get(2), vertices.get(3), null, 2.0), new Edge(
					vertices.get(3), vertices.get(0), null, 2.0), new Edge(vertices.get(
						1), vertices.get(3), null, 2.0 * Math.sqrt(2.0)), new Edge(vertices
							.get(2), vertices.get(4), null, 5.0), new Edge(vertices.get(0),
								vertices.get(5), null, 5.0), new Edge(vertices.get(4), vertices
									.get(6), null, 9.0), new Edge(vertices.get(5), vertices.get(
										6), null, 9.0));
		return createGraph(edges, vertices);
	}

	/**
	 * Creates a triangle graph with square cluster, three loops ("o"), two
	 * parallel edges and a dead end edge
	 *
	 * <pre>
	 *             o
	 *             4
	 *           _/||
	 *     3--2_/  ||
	 *     |\_|	   ||
	 *    _0--1	   ||
	 *  _/         ||
	 * /           ||
	 *5------------6------------7
	 *o------------o
	 * </pre>
	 */
	private static Graph createTriangleWithSquareClusterAndArtefacts() {
		final List<Vertex> vertices = Stream.generate(Vertex::new).limit(8).collect(
			Collectors.toList());
		vertices.get(0).addPoint(new Point(-1, -1, 0));
		vertices.get(1).addPoint(new Point(-1, 1, 0));
		vertices.get(2).addPoint(new Point(1, 1, 0));
		vertices.get(3).addPoint(new Point(1, -1, 0));
		vertices.get(4).addPoint(new Point(5, 4, 0));
		vertices.get(5).addPoint(new Point(-4, -5, 0));
		vertices.get(6).addPoint(new Point(5, -5, 0));
		vertices.get(7).addPoint(new Point(7, -5, 0));
		final List<Edge> edges = Arrays.asList(new Edge(vertices.get(0), vertices
			.get(1), null, 2.0), new Edge(vertices.get(1), vertices.get(2), null,
				2.0), new Edge(vertices.get(2), vertices.get(3), null, 2.0), new Edge(
					vertices.get(3), vertices.get(0), null, 2.0), new Edge(vertices.get(
						1), vertices.get(3), null, 2.0 * Math.sqrt(2.0)), new Edge(vertices
							.get(2), vertices.get(4), null, 5.0), new Edge(vertices.get(0),
								vertices.get(5), null, 5.0), new Edge(vertices.get(4), vertices
									.get(6), null, 9.0), new Edge(vertices.get(6), vertices.get(
										4), null, 9.0), // opposite-way-parallel-edge
			new Edge(vertices.get(5), vertices.get(6), null, 9.0), new Edge(vertices
				.get(5), vertices.get(6), null, 9.0), // same-way-parallel-edge
			new Edge(vertices.get(6), vertices.get(7), null, 2.0), // dead-end
			new Edge(vertices.get(5), vertices.get(5), null, 9.0), // loop
			new Edge(vertices.get(6), vertices.get(6), null, 9.0), // loop
			new Edge(vertices.get(4), vertices.get(4), null, 9.0) // loop
		);
		return createGraph(edges, vertices);
	}

	private static boolean hasPoint(final Edge edge, final Point point) {
		return edge.getV1().isVertexPoint(point) || edge.getV2().isVertexPoint(
			point);
	}

	private static boolean hasPointOnce(final Vertex vertex, final Point point) {
		return vertex.getPoints().stream().filter(p -> p.equals(point))
			.count() == 1;
	}
	// endregion
}
