
#include <mrc_json.h>
#include <json-builder.h>

#include <string.h>
#include <stdio.h>
#include <assert.h>

int
main(int argc, char **argv)
{
  char json_str[] = "{ \"i\": 2, \"d\": 2.3, \"i3\": [ 2, 3, 4] }";
  mrc_json_t json = mrc_json_parse(json_str);

  char *buf = mrc_json_to_string(json);
  printf("%s\n\n", buf);
  free(buf);
  
  mrc_json_print(json, 0);

  return 0;
}
