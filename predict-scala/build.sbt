name := "predict"
version := "0.1.0"
scalaVersion := "2.12.4"
organization := "com.hwaipy"
libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"
//libraryDependencies += "ch.obermuhlner" % "big-math" % "2.0.0" % "test"
scalacOptions ++= Seq("-feature")
scalacOptions ++= Seq("-deprecation")