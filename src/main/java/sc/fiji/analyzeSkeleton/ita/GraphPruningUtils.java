
package sc.fiji.analyzeSkeleton.ita;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.Vertex;

/**
 * Utility functions for pruning certain kinds of edges (and nodes) from a
 * {@link Graph}.
 *
 * @author Richard Domander
 * @author Alessandro Felder
 */
public final class GraphPruningUtils {

	private GraphPruningUtils() {}

	public static Graph pruneShortEdges(final Graph graph,
		final double minDistance, final boolean iterate, final boolean clustered)
	{
		return pruneShortEdges(graph, minDistance, iterate, clustered,
			new double[] { 1.0, 1.0, 1.0 });
	}

	public static Graph pruneShortEdges(final Graph graph,
		final double minDistance, final boolean iterate, final boolean clustered,
		final double[] voxelSize)
	{
		Graph pruned = graph.clone();
		boolean prune = true;
		removeLoops(pruned);
		pruned.getEdges().forEach(e -> euclideanDistance(e, voxelSize));
		while (prune) {
			final int startSize = pruned.getVertices().size();
			pruneDeadEnds(pruned, minDistance);
			if (clustered) {
				pruned = cleaningStep(pruned, minDistance, voxelSize);
			}
			else {
				pruned = edgeCleaning(pruned, minDistance, voxelSize);
			}
			removeParallelEdges(pruned);
			final int cleanedSize = pruned.getVertices().size();
			prune = iterate && startSize != cleanedSize;
		}
		return pruned;
	}

	/**
	 * Removes all loop edges from the graph.
	 *
	 * @param graph a graph.
	 * @see #isLoop(Edge)
	 */
	public static void removeLoops(final Graph graph) {
		final List<Edge> loops = graph.getEdges().stream().filter(
			GraphPruningUtils::isLoop).collect(toList());
		loops.forEach(GraphPruningUtils::removeBranchFromEndpoints);
		graph.getEdges().removeAll(loops);
	}

	/**
	 * Removes parallel edges from the graph, leaving at most one connection
	 * between each vertex pair.
	 * <p>
	 * While the input graph might already have parallel edges, creating them is a
	 * side effect how the graph is pruned in
	 * {@link #pruneDeadEnds(Graph, double)}. Note that the method is already
	 * called there, so there's no need for the user to explicitly call it
	 * afterwards.
	 * </p>
	 * <p>
	 * An edge is parallel, if there's another edge between its endpoint vertices.
	 * </p>
	 * <p>
	 * NB non-deterministic in choosing which of the parallel edges is kept.
	 * </p>
	 *
	 * @param graph A {@link Graph} that's assumed undirected
	 */
	public static void removeParallelEdges(final Graph graph) {
		final Map<Vertex, Integer> idMap = mapVertexIds(graph.getVertices());
		final Collection<Long> connections = new HashSet<>();
		final Collection<Edge> parallelEdges = new ArrayList<>();
		graph.getEdges().forEach(edge -> {
			final long hash = connectionHash(edge, idMap);
			if (!connections.add(hash)) {
				parallelEdges.add(edge);
			}
		});
		parallelEdges.forEach(GraphPruningUtils::removeBranchFromEndpoints);
		graph.getEdges().removeAll(parallelEdges);
	}

	// region -- Helper methods --

	private static double[] centroid(final Collection<Point> points) {
		final double[] centroid = new double[3];
		points.forEach(p -> {
			centroid[0] += p.x;
			centroid[1] += p.y;
			centroid[2] += p.z;
		});
		for (int i = 0; i < centroid.length; i++) {
			centroid[i] /= points.size();
		}
		return centroid;
	}

	private static Graph cleaningStep(final Graph graph, final double minDistance,
		final double[] voxelSize)
	{
		final List<Set<Vertex>> clusters = findClusters(graph, minDistance);
		final List<Vertex> clusterCentres = clusters.stream().map(
			GraphPruningUtils::getClusterCentre).collect(toList());
		final Map<Edge, Edge> replacements = mapReplacementEdges(clusters,
			clusterCentres);
		return createCleanGraph(graph, clusters, clusterCentres, replacements,
			voxelSize);
	}

	/**
	 * Creates a unique hash number for an edge based on its endpoints.
	 *
	 * @param e an edge in a graph.
	 * @param idMap mapping of vertices to unique ids.
	 * @return a hash code.
	 * @see #mapVertexIds(List)
	 */
	private static long connectionHash(final Edge e,
		final Map<Vertex, Integer> idMap)
	{
		final long nVertices = idMap.size();
		final long a = idMap.get(e.getV1());
		final long b = idMap.get(e.getV2());
		return a < b ? a * nVertices + b : b * nVertices + a;
	}

	private static Graph createCleanGraph(final Graph graph,
		final Collection<Set<Vertex>> clusters,
		final Collection<Vertex> clusterCentres, final Map<Edge, Edge> replacements,
		final double[] voxelSize)
	{
		final Collection<Edge> clusterEdges = getClusterEdges(replacements,
			voxelSize);
		final List<Edge> nonClusterEdges = graph.getEdges().stream().filter(
			e -> !replacements.containsKey(e) && isNotInClusters(e, clusters))
			.collect(toList());
		final Graph cleanGraph = new Graph();
		final Collection<Edge> cleanEdges = new HashSet<>();
		cleanEdges.addAll(nonClusterEdges);
		cleanEdges.addAll(clusterEdges);
		cleanGraph.getEdges().addAll(cleanEdges);
		clusterCentres.forEach(cleanGraph::addVertex);
		endpoints(nonClusterEdges).forEach(cleanGraph::addVertex);
		endpoints(clusterEdges).forEach(cleanGraph::addVertex);
		lonelyVertices(graph).forEach(cleanGraph::addVertex);
		removeDanglingBranches(cleanGraph);
		return cleanGraph;
	}

	private static Graph edgeCleaning(final Graph graph, final double minDistance,
		final double[] voxelSize)
	{
		Graph cleanGraph = graph.clone();
		final List<Edge> innerEdges = cleanGraph.getEdges().stream().filter(
			e -> isShort(e, minDistance) && !isDeadEnd(e)).collect(toList());
		for (final Edge innerEdge : innerEdges) {
			final List<Set<Vertex>> pairs = new ArrayList<>();
			pairs.add(getEndpoints(innerEdge));
			final List<Vertex> centroids = pairs.stream().map(
				GraphPruningUtils::getClusterCentre).collect(toList());
			final Map<Edge, Edge> replacements = mapReplacementEdges(pairs,
				centroids);
			final Graph tmp = createCleanGraph(cleanGraph, pairs, centroids,
				replacements, voxelSize);
			updateInnerEdges(innerEdges, replacements);
			cleanGraph = tmp;
		}
		return cleanGraph;
	}

	private static Stream<Vertex> endpoints(final Collection<Edge> edges) {
		return edges.stream().flatMap(e -> Stream.of(e.getV1(), e.getV2()))
			.distinct();
	}

	/**
	 * Sets the length of the {@link Edge} to the calibrated euclidean distance
	 * between its endpoints
	 */
	private static void euclideanDistance(final Edge e,
		final double[] voxelSize)
	{
		final double[] centre = centroid(e.getV1().getPoints());
		final double[] centre2 = centroid(e.getV2().getPoints());
		for (int i = 0; i < centre.length; i++) {
			centre[i] -= centre2[i];
		}
		final double l = length(centre, voxelSize);
		e.setLength(l);
	}

	/**
	 * Finds all the vertices in the cluster that has the given vertex.
	 * <p>
	 * A vertex is in the cluster if its connected to the start directly or
	 * indirectly via edges that have length less than the given distance.
	 * </p>
	 */
	private static Set<Vertex> fillCluster(final Vertex start,
		final double minDistance)
	{
		final Set<Vertex> cluster = new HashSet<>();
		final Stack<Vertex> stack = new Stack<>();
		stack.push(start);
		while (!stack.isEmpty()) {
			final Vertex vertex = stack.pop();
			cluster.add(vertex);
			final Set<Vertex> freeNeighbours = vertex.getBranches().stream().filter(
				e -> isShort(e, minDistance)).map(e -> e.getOppositeVertex(vertex))
				.filter(v -> !cluster.contains(v)).collect(toSet());
			stack.addAll(freeNeighbours);
		}
		return cluster;
	}

	/** Finds all the vertices that are in one of the graph's clusters */
	private static List<Vertex> findClusterVertices(final Graph graph,
		final double minDistance)
	{
		return graph.getEdges().stream().filter(e -> isShort(e, minDistance))
			.flatMap(e -> Stream.of(e.getV1(), e.getV2())).distinct().collect(
				toList());
	}

	private static Collection<Edge> getClusterEdges(
		final Map<Edge, Edge> replacements, final double[] voxelSize)
	{
		return replacements.values().stream().peek(e -> euclideanDistance(e,
			voxelSize)).collect(toList());
	}

	private static Set<Vertex> getEndpoints(final Edge edge) {
		final Set<Vertex> pair = new HashSet<>();
		pair.add(edge.getV1());
		pair.add(edge.getV2());
		return pair;
	}

	private static boolean isDeadEnd(final Edge e) {
		return Stream.of(e.getV1(), e.getV2()).filter(v -> v.getBranches()
			.size() == 1).count() == 1;
	}

	/**
	 * Checks if the edge forms a loop.
	 *
	 * @param edge an edge in a graph.
	 * @return true if both endpoints of the edge is the same vertex.
	 */
	private static boolean isLoop(final Edge edge) {
		return edge.getV1() != null && edge.getV1() == edge.getV2();
	}

	private static boolean isNotInClusters(final Edge e,
		final Collection<? extends Collection<Vertex>> clusters)
	{
		return clusters.stream().noneMatch(c -> c.contains(e.getV1()) && c.contains(
			e.getV2()));
	}

	private static double length(final double[] v, final double[] voxelSize) {
		final double x = v[0] * voxelSize[0];
		final double y = v[1] * voxelSize[1];
		final double z = v[2] * voxelSize[2];
		final double sqSum = DoubleStream.of(x, y, z).map(d -> d * d).sum();
		return Math.sqrt(sqSum);
	}

	private static Stream<Vertex> lonelyVertices(final Graph graph) {
		return graph.getVertices().stream().filter(v -> v.getBranches().isEmpty())
			.distinct();
	}

	/**
	 * Maps new replacement edges for the cluster centroid
	 * <p>
	 * When a cluster is pruned, that is, condensed to a single centroid vertex,
	 * all the edges that the cluster vertices had need to be connected to the new
	 * centroid. Instead of altering the edges, we create new ones to replace
	 * them, because {@link Edge} objects are immutable.
	 * <p>
	 * The method creates a map from old edges to their new replacements. Note
	 * that both endpoints of an edge in the mapping can change, when it's an edge
	 * connecting two clusters to each other.
	 * </p>
	 */
	private static Map<Edge, Edge> mapReplacementEdges(
		final List<Set<Vertex>> clusters, final List<Vertex> clusterCentres)
	{
		final Map<Edge, Edge> replacements = new HashMap<>();
		for (int i = 0; i < clusters.size(); i++) {
			final Collection<Vertex> cluster = clusters.get(i);
			final Vertex centre = clusterCentres.get(i);
			final Set<Edge> outerEdges = findEdgesWithOneEndInCluster(cluster);
			for (final Edge outerEdge : outerEdges) {
				final Edge oldEdge = replacements.getOrDefault(outerEdge, outerEdge);
				final Edge replacement = replaceEdge(oldEdge, cluster, centre);
				replacements.put(outerEdge, replacement);
			}
		}
		return replacements;
	}

	/**
	 * Assigns vertices a unique sequential id numbers.
	 *
	 * @param vertices list of vertices in a graph.
	 * @return a (vertex, id) mapping.
	 */
	private static Map<Vertex, Integer> mapVertexIds(
		final List<Vertex> vertices)
	{
		return IntStream.range(0, vertices.size()).boxed().collect(Collectors.toMap(
			vertices::get, Function.identity()));
	}

	private static void pruneDeadEnds(final Graph graph, final double tolerance) {
		final List<Edge> deadEnds = graph.getEdges().stream().filter(e -> isDeadEnd(
			e) && isShort(e, tolerance)).collect(toList());
		final List<Vertex> terminals = deadEnds.stream().flatMap(e -> Stream.of(e
			.getV1(), e.getV2())).filter(v -> v.getBranches().size() == 1).collect(
				toList());
		graph.getVertices().removeAll(terminals);
		deadEnds.forEach(GraphPruningUtils::removeBranchFromEndpoints);
		graph.getEdges().removeAll(deadEnds);
	}

	private static int[] realToIntegerCoordinate(final double[] centroid) {
		return Arrays.stream(centroid).mapToInt(d -> Double.isNaN(d)
			? Integer.MAX_VALUE : (int) Math.round(d)).toArray();
	}

	private static void removeBranchFromEndpoints(final Edge branch) {
		branch.getV1().getBranches().remove(branch);
		branch.getV2().getBranches().remove(branch);
	}

	/** Removes vertex branches that are no longer listed in the graph's edges */
	private static void removeDanglingBranches(final Graph graph) {
		graph.getVertices().stream().map(Vertex::getBranches).forEach(b -> b
			.removeIf(e -> !graph.getEdges().contains(e)));
	}

	private static Edge replaceEdge(final Edge edge,
		final Collection<Vertex> cluster, final Vertex centre)
	{
		final Vertex v1 = edge.getV1();
		final Vertex v2 = edge.getV2();
		final Edge replacement;
		if (cluster.contains(v1)) {
			replacement = new Edge(centre, v2, null, 0.0);
			replacement.getV1().setBranch(replacement);
			replacement.getV2().setBranch(replacement);
		}
		else if (cluster.contains(v2)) {
			replacement = new Edge(v1, centre, null, 0.0);
			replacement.getV1().setBranch(replacement);
			replacement.getV2().setBranch(replacement);
		}
		else {
			return null;
		}
		return replacement;
	}

	private static void updateInnerEdges(final List<Edge> innerEdges,
		final Map<Edge, Edge> replacements)
	{
		innerEdges.stream().filter(replacements::containsKey).forEach(e -> {
			final int i = innerEdges.indexOf(e);
			final Edge replacement = replacements.get(e);
			innerEdges.set(i, replacement);
		});
	}

	static List<Set<Vertex>> findClusters(final Graph graph,
		final double minDistance)
	{
		final List<Set<Vertex>> clusters = new ArrayList<>();
		final List<Vertex> clusterVertices = findClusterVertices(graph,
			minDistance);
		while (!clusterVertices.isEmpty()) {
			final Vertex start = clusterVertices.get(0);
			final Set<Vertex> cluster = fillCluster(start, minDistance);
			clusters.add(cluster);
			clusterVertices.removeAll(cluster);
		}
		return clusters;
	}

	/**
	 * Finds the edges that connect the cluster vertices to outside the cluster.
	 *
	 * @param cluster a collection of directly connected vertices.
	 * @return the edges that originate from the cluster but terminate outside it.
	 */
	static Set<Edge> findEdgesWithOneEndInCluster(
		final Collection<Vertex> cluster)
	{
		final Map<Edge, Long> edgeCounts = cluster.stream().flatMap(v -> v
			.getBranches().stream()).collect(groupingBy(Function.identity(),
				Collectors.counting()));
		return edgeCounts.keySet().stream().filter(e -> edgeCounts.get(e) == 1)
			.collect(toSet());
	}

	/**
	 * Creates a centroid vertex of all the vertices in a cluster.
	 *
	 * @param cluster a collection of directly connected vertices.
	 * @return A vertex at the geometric center of the cluster.
	 */
	static Vertex getClusterCentre(final Set<Vertex> cluster) {
		final Collection<Point> points = cluster.stream().flatMap(c -> c.getPoints()
			.stream()).collect(toList());
		final double[] centroid = centroid(points);
		final int[] coordinates = realToIntegerCoordinate(centroid);
		final Vertex vertex = new Vertex();
		vertex.addPoint(new Point(coordinates[0], coordinates[1], coordinates[2]));
		return vertex;
	}

	static boolean isShort(final Edge e, final double minLength) {
		return (e.getLength() < minLength);
	}
	// endregion
}
