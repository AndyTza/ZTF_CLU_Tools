### Tools for the Zwicky Transient Facility (ZTF) Census of the Local Universe (CLU) Experiment
___


![Aitoff Projection of Discovered ZTF CLU Transients](https://github.com/AndyTza/andytza.github.io/blob/master/images/CLU_Map.png?raw=true)


Before running the ZTF CLU scripts, you will need to add your fritz.science credentials. To start using the
fritz.scince API you will need to add to `/scripts/user_info.json` your last name and Command-Line Authentication Token:

```
{
  "user":
    {"user_last_name": "ADD_USER_LAST_NAME_HERE",
     "FritzToken": "ADD_FRITZ_TOKEN_HERE"
    }
}
```

#### Running the `ZTF CLU Luminous Annotation`
____
In ZTF Phase II the ZTF CLU experiment will be a volume-luminosity-limited supernovae sample in hope to find some of the most sub-luminous transients in the local universe.

This annotator will fetch the photometry of all sources saved to the Caltech ZTF CLU program (`program_id:43`) after the date indicated by the `-d` flag in format `YYYY-MM-DD`. The annotator will go through each source and annotate if they're luminous based on a specific luminosity cut (docs of `luminous_CLU_annotation.py` for more details). To run the ZTF CLU luminosity annotator, see the example below:

```
./luminous_CLU_annotation.py -d 2021-01-01
```

#### Running the `Fetch ZTF CLU Catalog`
___
In the case one would like to asses the general properties of the ZTF CLU experiment, the `fetch_CLU_souces.py` module is a quick method for downloading some basic parameters found in the ZTF CLU catalog (ascii table). The query returns the following parameters: `obj_id`, `saved_date`, `coordinates`, `redshift`, `CLU_d_gal` (seperation from CLU host in arcsec), `peak_absolute_magnitude`, `peak_apparent_magnitude`, `luminous_event`(boolean, if the source is luminous or sub-luminous give a -17 cut). To run the code, see example below:

```
./fetch_CLU_sources.py

```
`Note`: It roughly takes ~10 minutes to generate a table of ~3000 sources. If you would like more parameters to be added please make a PR or contact me. 
