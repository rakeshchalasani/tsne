package tsne

import java.io.File
import breeze.io.TextReader.FileReader
import breeze.linalg._
import breeze.numerics._
import breeze.plot._
import breeze.stats._
import breeze.stats.distributions._
import java.awt.Color
/**
 * Implementation of t-SNE ( t-distributed stochastic neighbourhood embedding).
 * ref: L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE.
 *      Journal of Machine Learning Research 9(Nov):2579-2605, 2008
 *
 */

class tSNE(val ndims: Int = 2, initDims:Int = 50) {

  private def Hbeta(D:DenseVector[Double], beta:Double = 1.0): (Double, DenseVector[Double]) = {
    // Compute the perplexity and the P-row for a specific value of the precision of a Gaussian distribution.

    // Compute P-row and corresponding perplexity
    val P = exp(- D * beta)
    val sumP = sum(P)
    val H = log(sumP) + (beta * sum(D :* P) / sumP)
    P := (P / sumP)
    (H, P)
  }

  private def x2p(X:DenseMatrix[Double], tol:Double, perplexity:Int) = {
    // Performs a binary search to get P-values in such a way that each conditional Gaussian has the same perplexity.

    // Initialize some variables
    print("Computing pairwise distances...")
    val n = X.cols
    val sum_X = { val t = pow(X,2); sum(t(::,*)).toDenseVector }
    // Pair-wise distance between each pair of points.
    val D = {
      val t = -2.0 :* (X.t * X)
      val x = (t(*,::) + sum_X)
      x(::,*) + sum_X
    }.asInstanceOf[DenseMatrix[Double]]
    val P = DenseMatrix.zeros[Double](n, n)
    val beta = DenseVector.ones[Double](n)
    val logU = log(perplexity)

    var i = 0
    while(i < n){
      // Print progress
      if((i+1) % 500 == 0) println("Computing P-values for point ", i+1, " of ", n, "...")

      // Compute the Gaussian kernel and entropy for the current precision
      var betamin = Double.NegativeInfinity
      var betamax =  Double.PositiveInfinity
      val t1 = D(0 until i,i)
      val t2 = D(i+1 until n, i)
      val Di = DenseVector.vertcat(t1, t2)
      var entP = Hbeta(Di, beta(i))

      // Evaluate whether the perplexity is within tolerance
      var Hdiff = entP._1 - logU
      var tries = 0
      while(abs(Hdiff) > tol && tries < 50) {
        // If not, increase or decrease precision
        if(Hdiff > 0) {
          if (betamax == Double.PositiveInfinity || betamax == Double.NegativeInfinity)
            beta(i) = beta(i) * 2
          else
            beta(i) = (beta(i) + betamax) / 2
          betamin = beta(i)
        } else {
          if(betamin == Double.PositiveInfinity || betamin == Double.NegativeInfinity)
            beta(i) = beta(i) / 2
          else
            beta(i) = (beta(i) + betamin) / 2
          betamax = beta(i)
        }

        // Recompute the values
        entP = Hbeta(Di, beta(i))
        Hdiff = entP._1 - logU
        tries = tries + 1
      }

      // Set the final row of P
      var j = 0
      while(j < n-1) {
        if (j < i) P(i, j) = entP._2(j)
        else if (j > i) P(i, j) = entP._2(j - 1)
        else Unit
        j += 1
      }

      i += 1
    }

    // Return final P-matrix
    println("Mean value of sigma: " +  mean(sqrt(1.0 / beta)))
    P
  }

  private def impl_tsne(data:DenseMatrix[Double], labels: DenseVector[Double], perplexity:Int):DenseMatrix[Double] = {
    val (d, n) = (data.rows, data.cols)
    // Initialize model.
    data := data - tile(mean(data(*,::)), 1, n)
    val X = PCA(initDims).compute(data)
    //    val X = data
    val max_iter = 1000
    val initial_momentum = 0.5
    val final_momentum = 0.8
    val eta:Double = 500
    val min_gain = 0.01
    val Y:DenseMatrix[Double] = DenseMatrix.rand(ndims, n, Gaussian(0,1))
    val dY = DenseMatrix.zeros[Double](ndims, n)
    val iY = DenseMatrix.zeros[Double](ndims, n)
    var gains = DenseMatrix.ones[Double](ndims, n)

    // Compute P-values
    val P = x2p(X, 1e-5, perplexity)
    P := P + P.t
    P := P :/ sum(P)
    P := P * 4.0									// early exaggeration
    P := max(P, 1e-12)

    val fig = Figure()
    val fillColors = Array(Color.RED, Color.BLACK, Color.BLUE, Color.MAGENTA, Color.YELLOW, Color.GRAY, Color.GREEN, Color.PINK, Color.DARK_GRAY, Color.CYAN)
    val pointColors = labels.map(n => fillColors(n.toInt))


    // Run iterations
    for(iter <- 0 to max_iter) {

      // Compute pairwise affinities
      val sum_Y = {val t = pow(Y,2); sum(t(::,*))}.toDenseVector
      val num = {
        val t = -2.0 * (Y.t * Y)
        val y = (t(*,::) + sum_Y)
        ((y(::,*) + sum_Y) + 1.0).map(x => 1.0 / x)
      }.asInstanceOf[DenseMatrix[Double]]

      //      num(0 until n, 0 until n) := 0.0

      var i = 0
      while(i < num.cols) {
        num(i, i) = 0.0
        i += 1
      }
      val Q:DenseMatrix[Double] = num / sum(num)
      Q := max(Q, 1e-12)

      // Compute gradient
      val PQ = P - Q
      for(i <- 0 until n) {
        val t1 = tile((PQ(i,::) :* num(i,::)).t,1, ndims)
        val t2 = (tile(Y(::, i), 1, n)-Y).t
        val t = sum(t1 :* t2, Axis._0).toDenseVector
        dY(::,i) := t
      }

      // Perform the update
      val momentum = if(iter < 20) initial_momentum else final_momentum

      val t = ((gains + 0.2) :* ((dY :> 0.0) :!= (iY :> 0.0)).map(x => if(x) 1.0 else 0.0))
      gains = t + ((gains * 0.8) :* ((dY :> 0.0) :== (iY :> 0.0)).map(x => if(x) 1.0 else 0.0))
      gains = gains.map(x => if(x < min_gain) min_gain else x)
      iY := (momentum :* iY) - (eta :* (gains :* dY))
      Y := Y + iY
      Y := Y - tile(mean(Y(*,::)),1,n)

      // Compute current value of cost function
      if((iter + 1) % 10 == 0) {
        val C = sum(P :* log(P :/ Q))
        println("Iteration " + (iter + 1) +  ": error is " +  C) //Stop lying about P -values
        fig.subplot(0) += scatter(Y(0,::).t, Y(1,::).t, (n:Int) => 1.0, colors = (n:Int) => pointColors(n))
      }
      if(iter == 100) P := P / 4.0
    }
    Y
  }

  /**
   * @param data (d x N) Dense matrix.
   * @return (ndims x N) Dense matrix.
   * */
  def compute(data:DenseMatrix[Double], labels: DenseVector[Double], perplexity:Int = 30):DenseMatrix[Double] = {
    impl_tsne(data, labels, perplexity)
  }

}

object TestTSNE extends App {

  //  val x = DenseMatrix.rand[Double](4,300, Rand.gaussian(0.0, 1.0))
  val x = ReadTSNE()
  val labels = ReadTSNELabels()
  val model = new tSNE(2, 50)

  model.compute(x, labels, 20)

}

case class PCA(ndims:Int) {

  def compute(data:DenseMatrix[Double]):DenseMatrix[Double] = {
    val xxt = data * data.t
    val pc = eig(xxt)
    pc.eigenvectors(::, 0 until ndims).t * data
  }
}

object ReadTSNE {
  def apply() = {

    val fileName = "/mnist/mnist2500_X.txt"
    val file = new File(getClass.getResource(fileName).toURI)

    val reader = new FileReader(file)
    var cont = true
    val line = reader.readRemaining().split("\n")
    val out = line.flatMap(x => {
      x.split(" ").filter(_ != "").map(_.toDouble)
    })
    new DenseMatrix[Double](28*28, 2500, out)
  }
}

object ReadTSNELabels {
  def apply() = {
    val fileName = "/mnist/mnist2500_labels.txt"
    val file = new File(getClass.getResource(fileName).toURI)

    val reader = new FileReader(file)
    val line = reader.readRemaining().split("\n")
    val out = line.flatMap(x => {
      x.split(" ").filter(_ != "").map(_.toDouble)
    })
    new DenseVector[Double](out)
  }
}