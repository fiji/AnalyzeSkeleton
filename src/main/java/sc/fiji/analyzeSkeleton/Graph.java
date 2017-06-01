/*
 * #%L
 * AnalyzeSkeleton_ plugin for ImageJ.
 * %%
 * Copyright (C) 2008 - 2017 Ignacio Arganda-Carreras.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package sc.fiji.analyzeSkeleton;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import java.util.ArrayList;
import java.util.Map;
import java.util.Stack;
import java.util.function.Function;

/**
 * This class represents an undirected graph to allow 
 * visiting the skeleton in an efficient way
 */
public class Graph 
{
	/** list of edges */
	private ArrayList < Edge > edges = null;
	/** list of vertices */
	private ArrayList < Vertex > vertices = null;

	/** root vertex */
	private Vertex root = null;

	// --------------------------------------------------------------------------
	/**
	 * Empty constructor.
	 */
	public Graph()
	{
		this.edges = new ArrayList<>();
		this.vertices = new ArrayList<>();
	}

	// --------------------------------------------------------------------------
	/**
	 * Add edge to the graph.
	 * @param e edge to be added
	 * @return false if the edge could not be added, true otherwise
	 */
	public boolean addEdge(Edge e)
	{
		if(this.edges.contains(e))
			return false;
		else
		{
			// Set vertices from e as neighbors (undirected graph)
			e.getV1().setBranch(e);
			if(! e.getV1().equals(e.getV2()))
				e.getV2().setBranch(e);
			// Add edge to the list of edges in the graph
			this.edges.add(e);
			return true;
		}
	}// end method addEdge

	// --------------------------------------------------------------------------
	/**
	 * Add vertex to the graph.
	 * @param v vertex to be added
	 * @return false if the vertex could not be added, true otherwise
	 */
	public boolean addVertex(Vertex v)
	{
		if(this.vertices.contains(v))
			return false;
		else
		{
			this.vertices.add(v);
			return true;
		}
	}// end method addVertex
	// --------------------------------------------------------------------------
	/**
	 * Get list of vertices in the graph.
	 * @return list of vertices in the graph
	 */
	public ArrayList<Vertex> getVertices()
	{
		return this.vertices;
	}
	// --------------------------------------------------------------------------
	/**
	 * Get list of edges in the graph.
	 * @return list of edges in the graph
	 */
	public ArrayList<Edge> getEdges()
	{
		return this.edges;
	}

	// --------------------------------------------------------------------------
	/**
	 * Set root vertex.
	 */
	void setRoot(Vertex v)
	{
		this.root = v;
	}
	// --------------------------------------------------------------------------
	/**
	 * Get root vertex.
	 * @return root vertex of the graph
	 */
	public Vertex getRoot()
	{
		return this.root;
	}

	// --------------------------------------------------------------------------
	/**
	 * Depth first search method to detect cycles in the graph.
	 * 
	 * @return list of BACK edges
	 */
	ArrayList<Edge> depthFirstSearch()
	{
		ArrayList<Edge> backEdges = new ArrayList<Edge>();

		// Create empty stack
		Stack<Vertex> stack = new Stack<Vertex>();

		// Mark all vertices as non-visited
		for(final Vertex v : this.vertices)
			v.setVisited(false);

		// Push the root into the stack
		stack.push(this.root);

		int visitOrder = 0;

		while(!stack.empty())
		{
			Vertex u = stack.pop();

			if(!u.isVisited())
			{
				//IJ.log(" Visiting vertex " + u.getPoints().get(0));

				// If the vertex has not been visited yet, then
				// the edge from the predecessor to this vertex
				// is mark as TREE
				if(u.getPredecessor() != null)
					u.getPredecessor().setType(Edge.TREE);

				// mark as visited
				u.setVisited(true, visitOrder++);

				for(final Edge e : u.getBranches())
				{
					// For the undefined branches:
					// We push the unvisited vertices in the stack,
					// and mark the edge to the others as BACK 
					if(e.getType() == Edge.UNDEFINED)
					{
						final Vertex ov = e.getOppositeVertex(u);
						if(!ov.isVisited())
						{
							stack.push(ov);
							ov.setPredecessor(e);
						}
						else
						{
							e.setType(Edge.BACK);
							backEdges.add(e);
						}

					}
				}
			}
		}

		return backEdges;

	} // end method depthFirstSearch

    @Override
	public Graph clone() {
		final Graph clone = new Graph();
		final Map<Vertex, Vertex> vertexMap = vertices.stream().collect(toMap(
			identity(), Vertex::cloneUnconnected));
		final Function<Edge, Edge> cloner = e -> e.clone(vertexMap.get(e.getV1()),
			vertexMap.get(e.getV2()));
		final Map<Edge, Edge> edgeMap = edges.stream().collect(toMap(identity(),
			cloner));
		vertices.forEach(v -> vertexMap.get(v).setPredecessor(edgeMap.get(v
			.getPredecessor())));
		// Iterate in the order of the originals to preserve order in the cloned
		// graph (makes testing easier)
		edges.stream().map(edgeMap::get).forEach(clone::addEdge);
		vertices.stream().map(vertexMap::get).forEach(clone::addVertex);
		clone.setRoot(vertexMap.get(root));
		return clone;
	}

}// end class Graph
