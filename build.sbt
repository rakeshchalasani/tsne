
// factor out common settings into a sequence
lazy val commonSettings = Seq(
  organization := "tsne",
  version := "0.1.0",
  // set the Scala version used for the project
  scalaVersion := "2.11.5"
)


lazy val root = (project in file(".")).
  settings(commonSettings: _*).
  settings(libraryDependencies ++= Seq(
      // other dependencies here
      "org.scalanlp" % "breeze_2.11" % "0.11-M0",
      // native libraries are not included by default. add this if you want them (as of 0.7)
      // native libraries greatly improve performance, but increase jar sizes.
      "org.scalanlp" % "breeze-natives_2.11" % "0.11-M0",
      "org.scalanlp" % "breeze-viz_2.11" % "0.11-M0"
    )
  )