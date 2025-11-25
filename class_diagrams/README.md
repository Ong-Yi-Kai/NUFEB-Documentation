Here we use [PlantUML](https://plantuml.com/) to build class diagrams. 

## Install PlantUML

### Linux
To install plantuml on Linux:
1. Install Java if not installed: `sudo apt install default-jre`
2. Install Graphviz: `sudo apt install graphviz`
3. Install PlantUML: `sudo apt install plantuml -y`



### Mac
To install plantuml on Mac
1. Make sure you have java: `java -version`
2. Install Graphviz: `brew install graphviz`
3. Install PlantUML: `brew install PlantUML`


## To create diagrams
In the command line run

```bash
plantuml -tpng diagram.puml  # for PNG
plantuml -tsvg diagram.puml  # for SVG
```
