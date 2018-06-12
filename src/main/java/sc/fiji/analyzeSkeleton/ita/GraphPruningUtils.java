
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
		final double tolerance)
	{
		Graph pruned = graph.clone();
		// final int prunedSize = pruned.getVertices().size();
		pruned.getEdges().forEach(GraphPruningUtils::euclideanDistance);
		pruneDeadEnds(pruned, tolerance);
		pruned = cleaningStep(pruned, tolerance);
		// final int cleanedSize = cleaned.getVertices().size();
		return pruned;
	}

	// region -- Helper methods =--
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

	private static Graph cleaningStep(final Graph graph, final double tolerance) {
		final List<Set<Vertex>> clusters = findClusters(graph, tolerance);
		final Map<Set<Vertex>, Vertex> clusterCentres = clusters.stream().collect(
			Collectors.toMap(Function.identity(),
				GraphPruningUtils::getClusterCentre));
		final Map<Edge, Edge> replacements = new HashMap<>();
		clusterCentres.forEach((cluster, centre) -> mapReplacementEdges(
			replacements, cluster, centre));
		final Collection<Edge> clusterConnectingEdges = replacements.values();
		clusterConnectingEdges.forEach(GraphPruningUtils::euclideanDistance);
		return createCleanGraph(graph, clusters, clusterCentres, replacements,
			clusterConnectingEdges);
	}

	private static Graph createCleanGraph(final Graph graph,
		final List<Set<Vertex>> clusters,
		final Map<Set<Vertex>, Vertex> clusterCentres,
		final Map<Edge, Edge> replacements,
		final Collection<Edge> clusterConnectingEdges)
	{
		final List<Edge> nonClusterEdges = graph.getEdges().stream().filter(
			e -> !replacements.containsKey(e) && isNotInClusters(e, clusters))
			.collect(toList());
		final Graph cleanGraph = new Graph();
		final Collection<Edge> cleanEdges = new HashSet<>();
		cleanEdges.addAll(nonClusterEdges);
		cleanEdges.addAll(clusterConnectingEdges);
		cleanGraph.getEdges().addAll(cleanEdges);
		clusterCentres.values().forEach(cleanGraph::addVertex);
		endpoints(nonClusterEdges).forEach(cleanGraph::addVertex);
		endpoints(clusterConnectingEdges).forEach(cleanGraph::addVertex);
		lonelyVertices(graph).forEach(cleanGraph::addVertex);
		removeDanglingBranches(cleanGraph);
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
	private static void euclideanDistance(final Edge e) {
		final double[] centre = centroid(e.getV1().getPoints());
		final double[] centre2 = centroid(e.getV2().getPoints());
		for (int i = 0; i < centre.length; i++) {
			centre[i] -= centre2[i];
		}
		final double l = length(centre);
		e.setLength(l);
	}

	/**
	 * Finds all the vertices in the cluster that has the given vertex
	 * <p>
	 * A vertex is in the cluster if its connected to the start directly or
	 * indirectly via edges that have length less than the given tolerance
	 * </p>
	 */
	private static Set<Vertex> fillCluster(final Vertex start,
		final Double tolerance)
	{
		final Set<Vertex> cluster = new HashSet<>();
		final Stack<Vertex> stack = new Stack<>();
		stack.push(start);
		while (!stack.isEmpty()) {
			final Vertex vertex = stack.pop();
			cluster.add(vertex);
			final Set<Vertex> freeNeighbours = vertex.getBranches().stream().filter(
				e -> isShort(e, tolerance)).map(e -> e.getOppositeVertex(vertex))
				.filter(v -> !cluster.contains(v)).collect(toSet());
			stack.addAll(freeNeighbours);
		}
		return cluster;
	}

	/** Finds all the vertices that are in one of the graph's clusters */
	private static List<Vertex> findClusterVertices(final Graph graph,
		final Double tolerance)
	{
		return graph.getEdges().stream().filter(e -> isShort(e, tolerance)).flatMap(
			e -> Stream.of(e.getV1(), e.getV2())).distinct().collect(toList());
	}

	private static boolean isDeadEnd(final Edge e) {
		return Stream.of(e.getV1(), e.getV2()).filter(v -> v.getBranches()
			.size() == 1).count() == 1;
	}

	private static boolean isNotInClusters(final Edge e,
		final Collection<Set<Vertex>> clusters)
	{
		return clusters.stream().noneMatch(c -> c.contains(e.getV1()) && c.contains(
			e.getV2()));
	}

	private static double length(final double[] v) {
		final double sqSum = Arrays.stream(v).map(d -> d * d).sum();
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

	private static void mapReplacementEdges(final Map<Edge, Edge> replacements,
		final Collection<Vertex> cluster, final Vertex centroid)
	{
		final Set<Edge> outerEdges = findEdgesWithOneEndInCluster(cluster);
		outerEdges.forEach(outerEdge -> {
			final Edge oldEdge = replacements.getOrDefault(outerEdge, outerEdge);
			final Edge replacement = replaceEdge(oldEdge, cluster, centroid);
			replacements.put(outerEdge, replacement);
		});
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

	static List<Set<Vertex>> findClusters(final Graph graph,
		final double tolerance)
	{
		final List<Set<Vertex>> clusters = new ArrayList<>();
		final List<Vertex> clusterVertices = findClusterVertices(graph, tolerance);
		while (!clusterVertices.isEmpty()) {
			final Vertex start = clusterVertices.get(0);
			final Set<Vertex> cluster = fillCluster(start, tolerance);
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
	static Vertex getClusterCentre(final Collection<Vertex> cluster) {
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
