package org.apache.spark.graphx

import scala.reflect.ClassTag
import scala.collection.mutable
import java.io._
import java.util.HashMap
//import it.unimi.dsi.fastutil.longs._

import scala.language.existentials

import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.graphx.PartitionStrategy._
import org.apache.spark.graphx.lib._
import org.apache.spark.graphx.impl._
import org.apache.spark.graphx.util._
import org.apache.spark.internal.Logging
import org.apache.spark.storage.StorageLevel
import org.apache.spark.util.collection.SortDataFormat
import org.apache.spark.rdd.RDD

object TriangleCounting extends Logging{

  def main(args: Array[String]): Unit = {
    println("############  Triangle Count #################")
    if (args.length < 2) {
      System.err.println(
        "Usage: TriangleCounting <file> <0 undirected or 1 directed> <No.ofEstimators> --numEPart=<num_edge_partitions> [other options]")
      System.exit(1)
    }

    val fname = args(0)
    val directed = args(1).toInt // directed or undirected
    val r = args(2).toInt // # of estimators
    val optionsList = args.drop(3).map { arg =>
      arg.dropWhile(_ == '-').split('=') match {
        case Array(opt, v) => (opt -> v)
        case _ => throw new IllegalArgumentException("Invalid argument: " + arg)
      }
    }

    val options = mutable.Map(optionsList: _*)

    val conf = new SparkConf()
    GraphXUtils.registerKryoClasses(conf)


    val numEPart = options.remove("numEPart").map(_.toInt).getOrElse {
      println("Set the number of edge partitions using --numEPart.")
      sys.exit(1)
    }

    val partitionStrategy: Option[PartitionStrategy] = options.remove("partStrategy")
    .map(PartitionStrategy.fromString(_))
    val edgeStorageLevel = options.remove("edgeStorageLevel")
    .map(StorageLevel.fromString(_)).getOrElse(StorageLevel.MEMORY_ONLY)
    val vertexStorageLevel = options.remove("vertexStorageLevel")
    .map(StorageLevel.fromString(_)).getOrElse(StorageLevel.MEMORY_ONLY)

    val sc = new SparkContext(conf.setAppName("TriangleCount(" + fname + ")"))

    val graph = GraphLoader.edgeListFile(sc, fname,
      canonicalOrientation = false,
      numEdgePartitions = numEPart,
      edgeStorageLevel = edgeStorageLevel,
      vertexStorageLevel = vertexStorageLevel)
          // TriangleCount requires the graph to be partitioned
          .partitionBy(partitionStrategy.getOrElse(RandomVertexCut))//.cache()

    val canonicalGraph = graph.mapEdges(e => 1).removeSelfEdges().convertToCanonicalEdges() 
    val n = graph.vertices.count().toInt // # of vertices
    val m = graph.edges.count().toInt // # of edges
    println("########### Edges "+ m + "; ######## vertices " + n)

    
    //val accum = sc.accumulator(0)//sum accum

    val startTime = System.currentTimeMillis
    val parTriangle = canonicalGraph.edges.mapPartitions{
           iterator => {
          //var deg = new Array[Int](n+1)//deg map for bulk vertices
          var degMap = new java.util.HashMap[VertexId, Int]()
          var ptable = new java.util.HashMap[(VertexId,Int), Array[Int]]()
          var ltable = new java.util.HashMap[Edge[Int], Array[Int]]()
          var qtable = new java.util.HashMap[(VertexId,VertexId), Array[Int]]()
          var r1 = new Array[Edge[Int]](r)
          var r2 = new Array[Edge[Int]](r)
          var r2locs = new Array[Int](r)
          var betax = new Array[Int](r)
          var betay = new Array[Int](r)
          var c = new Array[Int](r)
          var close = new Array[Boolean](r)
          var locs = new Array[VertexId](r)


          val edgeList = iterator.toArray
          val w = edgeList.size
          val ran = scala.util.Random

          //Count type1 triangles for directed graph
          // if(directed != 0){
            

          // }
          // Step 1
          for( i <- 0 to (r-1)) {
            val loc = ran.nextInt(w)
            r1(i) = edgeList(loc)
            //c(i) = 0
            val tmpbuf = ltable.get(r1(i))
            if(tmpbuf == null){
              ltable.put(r1(i),Array(i))
              }else{
                ltable.put(r1(i),tmpbuf :+ i)
              }
            }

            edgeList.map(edge => {
              val src = edge.srcId
              val dst = edge.dstId

              val srcVal = degMap.get(src)
              val dstVal = degMap.get(dst)

              if(srcVal == null){
                degMap.put(src,1)
                }else
                {
                  degMap.put(src,srcVal+1)
                }

              if(directed == 0){
                if(dstVal == null){
                  degMap.put(dst,1)
                }else
                  {
                  degMap.put(dst,dstVal+1)
                  }
              }

                  val tmpbuf = ltable.get(edge)
                  if(tmpbuf!=null){
                    for( estimator <- tmpbuf) {
                      if(directed == 0){
                        betax(estimator) = degMap.get(src)
                      }
                      betay(estimator) = degMap.get(dst)
                    }           
                  }
                  })

          //step 2
          for( i <- 0 to (r-1)) {
            if(r1(i)!=null){
              val src = r1(i).srcId
              val dst = r1(i).dstId
              val a = degMap.get(src) - betax(i)
              val b = degMap.get(dst) - betay(i)
              if(directed == 0){
                c(i) = (a + b)
                if(c(i)>0){
                val alpha = 1+ran.nextInt(c(i))

                if(alpha <= a){
                  val tmpDeg = betax(i)+alpha
                  val tmpbuf = ptable.get((src,tmpDeg))
                  if(tmpbuf == null){
                    ptable.put((src,tmpDeg),Array(i))
                    }else{
                      ptable.put((src,tmpDeg),tmpbuf :+ i)
                    }
                    }else{
                      val tmpDeg = betax(i)+alpha-a
                      val tmpbuf = ptable.get((dst,tmpDeg))
                      if(tmpbuf == null){
                        ptable.put((dst,tmpDeg),Array(i))
                        }else{
                          ptable.put((dst,tmpDeg),tmpbuf :+ i)
                        }
                      }              
                    }
              }else{
                c(i) = b
                if(c(i)>0){
                  val alpha = 1+ran.nextInt(c(i))
                  val tmpDeg = betay(i)+alpha
                  val tmpbuf = ptable.get((dst,tmpDeg))
                  if(tmpbuf == null){
                      ptable.put((dst,tmpDeg),Array(i))
                        }else{
                      ptable.put((dst,tmpDeg),tmpbuf :+ i)
                        }
                }
              }
                  } 
                }

                degMap.clear()

                for( k <- 0 to (w-1) ) {
                  val src = edgeList(k).srcId
                  val dst = edgeList(k).dstId
                  val srcVal = degMap.get(src)
                  val dstVal = degMap.get(dst)

                  if(srcVal == null){
                    degMap.put(src,1)
                    }else
                    {
                      degMap.put(src,srcVal+1)
                    }
                    if(directed == 0){
                      if(dstVal == null){
                        degMap.put(dst,1)
                        }else
                        {
                          degMap.put(dst,dstVal+1)
                        }
                      }


                      val tmpbufx = ptable.get((src,degMap.get(src)))
                      val tmpbufy = ptable.get((dst,degMap.get(dst)))
                      if(directed == 0){
                        if(tmpbufx != null){

                        for(estimator <- tmpbufx) {

                          r2(estimator) = edgeList(k)
                          locs(estimator) = src
                          r2locs(estimator) = k
                        }
                        }else if(tmpbufy != null){

                          for(estimator <- tmpbufy) {

                            r2(estimator) = edgeList(k)
                            locs(estimator) = dst
                            r2locs(estimator) = k
                          }
                        }
                      }else{
                        if(tmpbufx != null){
                          for(estimator <- tmpbufx) {

                          r2(estimator) = edgeList(k)
                          //locs(estimator) = src
                          r2locs(estimator) = k
                        }
                        }

                      }


                      }

            //step 3
            for( i <- 0 to (r-1)) {
              var node1 = (-1).toLong
              var node2 = (-1).toLong
              
              if(directed == 0){
                if(r1(i) != null){
                  if(r1(i).srcId == locs(i))
                  node1 = r1(i).dstId
                  else
                  node1 = r1(i).srcId 
                }
                if(r2(i) != null){
                  if(r2(i).srcId == locs(i))
                  node2 = r2(i).dstId
                  else
                  node2 = r2(i).srcId
                }
              }else{
                if(r2(i) != null){
                node1 = r2(i).dstId
              }
               if(r1(i) != null){
                node2 = r1(i).srcId
              }
              }

              
              val tmpbuf = qtable.get((node1,node2))
              if(tmpbuf == null){
                qtable.put((node1,node2),Array(i))
                }else{
                  qtable.put((node1,node2),tmpbuf :+ i)
                }
              }

              for( k <- 0 to (w-1)) {
                val src = edgeList(k).srcId
                val dst = edgeList(k).dstId
                val tmpbuf = qtable.get((src,dst))

                if(directed ==0 ){
                  if(tmpbuf == null){
                    val tmpbuf2 = qtable.get((dst,src))
                    if(tmpbuf2 !=null){
                      for( estimator <- tmpbuf2) {
                       if(k>r2locs(estimator))
                       close(estimator) = true
                     } 
                   }  
                   }else{
                    // val splitStr = tmpbuf.trim.split("\\s+")
                    for( estimator <- tmpbuf) {
                      if(k>r2locs(estimator))
                      close(estimator) = true    
                    }

                  } 
                }else{
                  if(tmpbuf != null){
                     for( estimator <- tmpbuf) {
                       if(k>r2locs(estimator))
                       close(estimator) = true
                     } 
                  }

                }

                  }

                  var sum = 0
                  for( i <- 0 to (r-1)) {
                    if(close(i))
                    sum += c(i)  
                  }
                  val totalT = sum*w.toDouble/r.toDouble
                  Iterator(totalT)
                }
              }

              val totalTriangle = parTriangle.reduce((a, b) => a + b) * numEPart.toDouble * numEPart.toDouble 
              val totalTime = System.currentTimeMillis - startTime
              println("############# The number of triangles is " + totalTriangle)
              println("############# Total running time is " + totalTime.toDouble/1000.0)
              val pw = new PrintWriter(new File("output-"+System.currentTimeMillis.toDouble/1000.0/3600.0+".txt" ))
              pw.write("Triangles:"+totalTriangle+" Time:"+totalTime.toDouble/1000.0)
              pw.close
              sc.stop()

            }
          }
// scalastyle:on println
