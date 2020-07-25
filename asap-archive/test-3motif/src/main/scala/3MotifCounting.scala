package org.apache.spark.graphx

import scala.reflect.ClassTag
import scala.collection.mutable
import java.io._
//import java.util.HashMap
import it.unimi.dsi.fastutil.longs._
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap
import it.unimi.dsi.fastutil.ints._

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

object MotifCounting extends Logging{

  def main(args: Array[String]): Unit = {
    println("############  Motif Count #################")
    if (args.length < 1) {
      System.err.println(
        "Usage: MotifCounting <file> <No.ofEstimators> --numEPart=<num_edge_partitions> [other options]")
      System.exit(1)
    }

    val fname = args(0)
    //val fname = "s3n://graphx-datasets/twitter/part-0000[1-2]"
    val r = args(1).toInt // # of estimators
    val optionsList = args.drop(2).map { arg =>
      arg.dropWhile(_ == '-').split('=') match {
        case Array(opt, v) => (opt -> v)
        case _ => throw new IllegalArgumentException("Invalid argument: " + arg)
      }
    }

    val options = mutable.Map(optionsList: _*)

    // val conf = new SparkConf()
    val conf = new SparkConf().setMaster("local[1]")
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

    val sc = new SparkContext(conf.setAppName("MotifCount(" + fname + ")"))
    val graph = GraphLoader.edgeListFile(sc, fname,
      canonicalOrientation = true,
      numEdgePartitions = numEPart,
      edgeStorageLevel = edgeStorageLevel,
      vertexStorageLevel = vertexStorageLevel)
          // TriangleCount requires the graph to be partitioned
          .partitionBy(partitionStrategy.getOrElse(RandomVertexCut)).cache()

    // val n = graph.vertices.count().toInt // # of vertices
    // val m = graph.edges.count().toInt // # of edges
    // println("########### Edges "+ m + "; ######## vertices " + n)


    val startTime = System.currentTimeMillis
    val parTriangle = graph.edges.mapPartitions{
      iterator => {
          //var deg = new Array[Int](n+1)//deg map for bulk vertices
          val degMap = new Long2IntOpenHashMap()//fastutil hashmaps
          val ptable = new Object2ObjectOpenHashMap[(VertexId,Int), Array[Int]]()
          val ltable = new Long2ObjectOpenHashMap[Array[Int]]()
          //val qtable = new Object2ObjectOpenHashMap[(VertexId,VertexId), Array[Int]]()

          val r1 = new Array[Edge[Int]](r)
          //val r2 = new Array[Edge[Int]](r)
          //val r2locs = new Array[Int](r)

          val betax = new Array[Int](r)
          val betay = new Array[Int](r)

          //val betaa = new Array[Int](r)
          //val betab = new Array[Int](r)

          val c = new Array[Int](r)
          //val close = new Array[Boolean](r)
          
          //val node1s = new Array[VertexId](r)
          //val node2s = new Array[VertexId](r)


          val edgeList = iterator.toArray
          val w = edgeList.size
          val ran = scala.util.Random

          // Step 1 sample first edge
          for( i <- 0 to (r-1)) {
            val loc = ran.nextInt(w)
            r1(i) = edgeList(loc)
            val tmpbuf = ltable.get(loc)
            if(tmpbuf == null){
              ltable.put(loc,Array(i))
              }else{
                ltable.put(loc,tmpbuf :+ i)
              }
            }

            for( k <- 0 to (w-1) ) {
              val src = edgeList(k).srcId
              val dst = edgeList(k).dstId

              val srcVal = degMap.get(src)
              val dstVal = degMap.get(dst)
              if(srcVal == null){
                degMap.put(src,1)
                }else
                {
                  degMap.addTo(src,1)
                }
              if(dstVal == null){
                  degMap.put(dst,1)
                }else
                  {
                    degMap.addTo(dst,1)
                  }

                  val tmpbuf = ltable.get(k)
                  if(tmpbuf!=null){
                    for( estimator <- tmpbuf) {
                      betax(estimator) = degMap.get(src)
                      betay(estimator) = degMap.get(dst)
                    }           
                  }
                  }
          //step 2 sample second edge
          for( i <- 0 to (r-1)) {
            if(r1(i)!=null){
              val src = r1(i).srcId
              val dst = r1(i).dstId
              val a = degMap.get(src) - betax(i)
              val b = degMap.get(dst) - betay(i)
              if((a+b)>0){        
                c(i) = (a + b)
                }else{
                  c(i) = 0
                }
              } 
            }

                //degMap.clear()

                  var sum = 0.0
                  for( i <- 0 to (r-1)) {
                    sum += c(i)
                  }
                  val totalT = sum*w.toDouble/r.toDouble
                  Iterator(totalT)
                }
              }

              val totalTriangle = parTriangle.reduce((a, b) => a + b) * numEPart.toDouble//add up in reduce phase
              val totalTime = System.currentTimeMillis - startTime
              println("############# The number of 3-Motifs is " + totalTriangle)
              println("############# Total running time is " + totalTime.toDouble/1000.0)

              val pw = new PrintWriter(new FileOutputStream(new File("output3Motif.txt"),true))
              pw.println(r+"  "+totalTriangle+"  "+totalTime.toDouble/1000.0+" "+numEPart)
              pw.close

              sc.stop()

            }
          }
// scalastyle:on println
