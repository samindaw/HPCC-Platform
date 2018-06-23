<?xml version="1.0" encoding="UTF-8"?>
<!--
##############################################################################
# HPCC SYSTEMS software Copyright (C) 2012 HPCC Systems®.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################
-->
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsd="http://www.w3.org/2001/XMLSchema">
    <xsl:output method="text" version="1.0" encoding="UTF-8" indent="yes"/>
    <xsl:param name="sourceFileName" select="'UNKNOWN'"/>
    <xsl:param name="importsList" select="''"/>
    <xsl:variable name="docname" select="/expesdl/esxdl/@name"/>
    <xsl:template match="/">
        <xsl:apply-templates select="expesdl"/>
    </xsl:template>
    <xsl:template match="text()" name="outputImports">
    <xsl:param name="import" select="."/>
        <xsl:if test="string-length($import) > 0">
            <xsl:text disable-output-escaping="yes">IMPORT </xsl:text><xsl:value-of select="$import"/><xsl:text disable-output-escaping="yes">;
</xsl:text>
        </xsl:if>
    </xsl:template>
    <xsl:template match="text()" name="handleImports">
        <xsl:param name="imports" select="."/>

        <xsl:choose>
            <xsl:when test="contains($imports, ',')">
                <xsl:call-template name="outputImports">
                    <xsl:with-param name="import" select="substring-before($imports, ',')"/>
                </xsl:call-template>
                <xsl:call-template name="handleImports">
                    <xsl:with-param name="imports" select="substring-after($imports, ',')"/>
                </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
                <xsl:call-template name="outputImports">
                    <xsl:with-param name="import" select="$imports"/>
                </xsl:call-template>
                <xsl:text disable-output-escaping="yes">
</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

    <xsl:template name="doNotChangeManuallyComment">
        <xsl:text disable-output-escaping="yes">/*** Not to be hand edited (changes will be lost on re-generation) ***/
/*** ECL Interface generated by esdl2ecl version 1.0 from </xsl:text><xsl:copy-of select="$sourceFileName"/> <xsl:text disable-output-escaping="yes">.xml. ***/
/*===================================================*/

</xsl:text>
    </xsl:template>
    <xsl:template match="expesdl">
            <xsl:apply-templates select="esxdl"/>
    </xsl:template>
    <xsl:template match="esxdl">
                 <xsl:call-template name="doNotChangeManuallyComment"/>
                  <xsl:call-template name="handleImports">
                    <xsl:with-param name="imports" select="$importsList"/>
                  </xsl:call-template>
            <xsl:text disable-output-escaping="yes">EXPORT </xsl:text>
			<xsl:choose>
				<xsl:when test="starts-with(@name, 'wsm_')"><xsl:value-of select="substring(@name, 5)"/></xsl:when>
				<xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
			</xsl:choose>
			<xsl:text disable-output-escaping="yes"> := MODULE

</xsl:text>
		<xsl:apply-templates select="EsdlStruct"/>
		<xsl:apply-templates select="EsdlRequest"/>
		<xsl:apply-templates select="EsdlResponse"/>

		<xsl:text disable-output-escaping="yes">
END;

</xsl:text>
          <xsl:call-template name="doNotChangeManuallyComment"/>
	</xsl:template>

	<xsl:template match="EsdlElement[@complex_type]">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
		  <xsl:text disable-output-escaping="yes">	</xsl:text><xsl:call-template name="output_ecl_complex_type"/><xsl:text disable-output-escaping="yes"> </xsl:text> <xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> {</xsl:text><xsl:call-template name="output_xpath"/><xsl:text disable-output-escaping="yes">};</xsl:text><xsl:call-template name="output_comments"/>
          <xsl:text disable-output-escaping="yes">
</xsl:text>
		</xsl:if>
	</xsl:template>
	<xsl:template match="EsdlElement[@type]">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
           <xsl:text disable-output-escaping="yes">	</xsl:text><xsl:call-template name="output_basic_type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> {</xsl:text><xsl:call-template name="output_xpath"/><xsl:if test="@ecl_max_len"><xsl:text disable-output-escaping="yes">, MAXLENGTH(</xsl:text><xsl:value-of select="@ecl_max_len"/><xsl:text disable-output-escaping="yes">)</xsl:text></xsl:if><xsl:text disable-output-escaping="yes">};</xsl:text><xsl:call-template name="output_comments"/>
<xsl:text disable-output-escaping="yes">
</xsl:text>
		</xsl:if>
	</xsl:template>
	<xsl:template match="EsdlEnum">
		<xsl:variable name="entype" select="@enum_type"/>
		<xsl:text disable-output-escaping="yes">	</xsl:text><xsl:call-template name="output_enum_type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> {</xsl:text><xsl:call-template name="output_xpath"/><xsl:text disable-output-escaping="yes">}; </xsl:text><xsl:call-template name="output_comments"/>
<xsl:text disable-output-escaping="yes">
</xsl:text>
	</xsl:template>
	<xsl:template match="EsdlArray[@type='string']|EsdlList[@type='string']">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
        <xsl:text disable-output-escaping="yes">	SET OF STRING </xsl:text><xsl:call-template name="output_ecl_name"/>
        <xsl:text disable-output-escaping="yes"> {XPATH('</xsl:text>
		<xsl:choose>
			<xsl:when test="name(.) ='EsdlArray'">
				<xsl:if test="not(@flat_array)"><xsl:value-of select="@name"/></xsl:if><xsl:text disable-output-escaping="yes">/</xsl:text><xsl:call-template name="output_item_tag"/>
			</xsl:when>
			<xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
		</xsl:choose>
                <xsl:text disable-output-escaping="yes">')</xsl:text>
		<xsl:choose>
			<xsl:when test="@max_count_var"><xsl:text disable-output-escaping="yes">, MAXCOUNT(</xsl:text><xsl:value-of select="@max_count_var"/><xsl:text disable-output-escaping="yes">)};</xsl:text></xsl:when>
			<xsl:when test="@max_count"><xsl:text disable-output-escaping="yes">, MAXCOUNT(</xsl:text><xsl:value-of select="@max_count"/><xsl:text disable-output-escaping="yes">)};</xsl:text></xsl:when>
		    <xsl:otherwise>
                <xsl:text disable-output-escaping="yes">, MAXCOUNT(1)}; // max_count must be specified in ESDL defintion! </xsl:text>
                <xsl:message terminate="no">EsdlArray MUST SPECIFY max_count</xsl:message>
		    </xsl:otherwise>
		</xsl:choose>
		<xsl:call-template name="output_comments"/>
		<xsl:text disable-output-escaping="yes">&#xa;</xsl:text>
	</xsl:if>
	</xsl:template>
	<xsl:template match="EsdlArray[starts-with(@ecl_type,'string') or starts-with(@ecl_type,'unicode')]">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
          <xsl:text disable-output-escaping="yes">&#9;</xsl:text>
		  <xsl:value-of select="@ecl_type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:call-template name="output_ecl_name"/>
          <xsl:text disable-output-escaping="yes"> {XPATH('</xsl:text>
		  <xsl:value-of select="@name"/>
		  <xsl:text disable-output-escaping="yes">')};</xsl:text>
		  <xsl:call-template name="output_comments"/>
		  <xsl:text disable-output-escaping="yes">&#xa;</xsl:text>
        </xsl:if>
	</xsl:template>
	<xsl:template match="EsdlArray|EsdlList">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
            <xsl:text disable-output-escaping="yes">	DATASET(</xsl:text> <xsl:call-template name="output_ecl_array_type"/><xsl:text disable-output-escaping="yes">) </xsl:text><xsl:call-template name="output_ecl_name"/>
                    <xsl:text disable-output-escaping="yes"> {XPATH('</xsl:text>
			<xsl:choose>
				<xsl:when test="name(.) ='EsdlArray'">
 					<xsl:if test="not(@flat_array)"><xsl:value-of select="@name"/></xsl:if>
					<xsl:text disable-output-escaping="yes">/</xsl:text>
					<xsl:call-template name="output_item_tag"/>
                                </xsl:when>
                        	<xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
			</xsl:choose>
            		<xsl:text disable-output-escaping="yes">')</xsl:text>
			<xsl:choose>
				<xsl:when test="@max_count_var"><xsl:text disable-output-escaping="yes">, MAXCOUNT(</xsl:text><xsl:value-of select="@max_count_var"/><xsl:text disable-output-escaping="yes">)};</xsl:text></xsl:when>
				<xsl:when test="@max_count"><xsl:text disable-output-escaping="yes">, MAXCOUNT(</xsl:text><xsl:value-of select="@max_count"/><xsl:text disable-output-escaping="yes">)};</xsl:text></xsl:when>
				<xsl:otherwise>
					<xsl:text disable-output-escaping="yes">, MAXCOUNT(1)}; // max_count must be specified in ESDL defintion! </xsl:text>
					<xsl:message terminate="no">EsdlArray MUST SPECIFY max_count</xsl:message>
				</xsl:otherwise>
			</xsl:choose>
		<xsl:call-template name="output_comments"/>
		<xsl:text disable-output-escaping="yes">
</xsl:text>
		</xsl:if>
	</xsl:template>

	<xsl:template match="EsdlStruct">
		<xsl:if test="not(@ecl_hide) and (@ecl_keep or not(@get_data_from))">
            <xsl:text disable-output-escaping="yes">EXPORT t_</xsl:text><xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> := RECORD</xsl:text><xsl:if test="@RecordCode"><xsl:text disable-output-escaping="yes"> //RecordCode[</xsl:text><xsl:value-of select="@RecordCode"/><xsl:text disable-output-escaping="yes">]</xsl:text></xsl:if>
			<xsl:if test="@base_type"><xsl:call-template name="output_ecl_base_type"/></xsl:if>
			<xsl:if test="@max_len"><xsl:text disable-output-escaping="yes">, MAXLENGTH (</xsl:text><xsl:value-of select="@max_len"/><xsl:text disable-output-escaping="yes">)</xsl:text></xsl:if>
			<xsl:text disable-output-escaping="yes">
</xsl:text>
			<xsl:apply-templates select="*"/>
            <xsl:if test="@element and not(*[@name='Content_'])">	STRING Content_ {XPATH('')};
</xsl:if>
            <xsl:text disable-output-escaping="yes">END;

</xsl:text>
		</xsl:if>
	</xsl:template>

	<xsl:template match="EsdlRequest">
        <xsl:if test="not(node())">
        <xsl:text disable-output-escaping="yes">/*Empty record generated from empty EsdlRequest
</xsl:text>
        </xsl:if>
        <xsl:text disable-output-escaping="yes">EXPORT t_</xsl:text><xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> := RECORD</xsl:text>
		<xsl:if test="@base_type"><xsl:call-template name="output_ecl_base_type"/></xsl:if>
		<xsl:if test="@max_len"><xsl:text disable-output-escaping="yes">, MAXLENGTH (</xsl:text><xsl:value-of select="@max_len"/><xsl:text disable-output-escaping="yes">)</xsl:text></xsl:if>
		<xsl:text disable-output-escaping="yes">
</xsl:text>
		<xsl:apply-templates select="*"/>
        <xsl:text disable-output-escaping="yes">END;
</xsl:text>
        <xsl:if test="not(node())">
        <xsl:text disable-output-escaping="yes">*/
</xsl:text>
        </xsl:if>
<xsl:text disable-output-escaping="yes">
</xsl:text>
	</xsl:template>

	<xsl:template match="EsdlResponse">
        <xsl:if test="not(node())">
        <xsl:text disable-output-escaping="yes">/*Empty record generated from empty EsdlResponse
</xsl:text>
        </xsl:if>
        <xsl:text disable-output-escaping="yes">EXPORT t_</xsl:text><xsl:call-template name="output_ecl_name"/><xsl:text disable-output-escaping="yes"> := RECORD</xsl:text>
		<xsl:if test="@base_type"><xsl:call-template name="output_ecl_base_type"/></xsl:if>
		<xsl:if test="@max_len"><xsl:text disable-output-escaping="yes">, MAXLENGTH (</xsl:text><xsl:value-of select="@max_len"/><xsl:text disable-output-escaping="yes">)</xsl:text></xsl:if>
		<xsl:text disable-output-escaping="yes">
</xsl:text>
		<xsl:apply-templates select="*"/>
        <xsl:text disable-output-escaping="yes">END;
</xsl:text>
        <xsl:if test="not(node())">
        <xsl:text disable-output-escaping="yes">*/
</xsl:text>
        </xsl:if>
<xsl:text disable-output-escaping="yes">
</xsl:text>
	</xsl:template>


<xsl:template name="output_ecl_name">

<xsl:choose>
	<xsl:when test="@ecl_name"><xsl:value-of select="@ecl_name"/></xsl:when>
	<xsl:otherwise>
		<xsl:variable name="nameword" select="translate(@name, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>
		<xsl:if test="/expesdl/keywords/keyword[@word=$nameword]"><xsl:text disable-output-escaping="yes">_</xsl:text></xsl:if><xsl:value-of select="@name"/>
	</xsl:otherwise>
</xsl:choose>

</xsl:template>

<xsl:template name="output_basic_type">
	<xsl:param name="basic_type" select="@type"/>
	<xsl:param name="size" select="@max_len"/>
	<xsl:choose>
		<xsl:when test="@ecl_type"><xsl:value-of select="@ecl_type"/><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:when>
        <xsl:when test="$basic_type='int'"><xsl:text disable-output-escaping="yes">INTEGER</xsl:text><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:when>
        <xsl:when test="$basic_type='long'"><xsl:text disable-output-escaping="yes">INTEGER4</xsl:text><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:when>
        <xsl:when test="$basic_type='short'"><xsl:text disable-output-escaping="yes">INTEGER2</xsl:text></xsl:when>
        <xsl:when test="$basic_type='int64'"><xsl:text disable-output-escaping="yes">INTEGER8</xsl:text></xsl:when>
        <xsl:when test="$basic_type='bool'"><xsl:text disable-output-escaping="yes">BOOLEAN</xsl:text></xsl:when>
        <xsl:when test="$basic_type='string'"><xsl:text disable-output-escaping="yes">STRING</xsl:text><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:when>
        <xsl:when test="$basic_type='double'"><xsl:text disable-output-escaping="yes">REAL8</xsl:text></xsl:when>
        <xsl:when test="$basic_type='float'"><xsl:text disable-output-escaping="yes">REAL4</xsl:text></xsl:when>
        <xsl:when test="$basic_type='base64Binary'"><xsl:text disable-output-escaping="yes">STRING</xsl:text></xsl:when>
		<xsl:when test="$basic_type"><xsl:value-of select="$basic_type"/><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:when>
        <xsl:otherwise><xsl:text disable-output-escaping="yes">STRING</xsl:text><xsl:if test="not(@ecl_max_len)"><xsl:value-of select="$size"/></xsl:if></xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="output_enum_type">
	<xsl:variable name="etype" select="@enum_type"/>
	<xsl:choose>
		<xsl:when test="/expesdl/types/type[@name=$etype]/@base_type">
			<xsl:call-template name="output_basic_type">
				<xsl:with-param name="basic_type" select="/expesdl/types/type[@name=$etype]/@base_type"/>
			</xsl:call-template>
		</xsl:when>
		<xsl:otherwise>
			<xsl:call-template name="output_basic_type"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>



<xsl:template name="output_ecl_base_type">
	<xsl:variable name="btype" select="@base_type"/>
	<xsl:if test="@base_type">
		<xsl:variable name="srcfile" select="/expesdl/types/type[@name=$btype]/@src"/>
		<xsl:text disable-output-escaping="yes"> (</xsl:text>
		<xsl:if test="$sourceFileName != $srcfile">
			<xsl:value-of select="$srcfile"/><xsl:text disable-output-escaping="yes">.</xsl:text>
		</xsl:if>
		<xsl:text disable-output-escaping="yes">t_</xsl:text><xsl:value-of select="$btype"/><xsl:text disable-output-escaping="yes">)</xsl:text>
	</xsl:if>
</xsl:template>

<xsl:template name="output_ecl_complex_type">
	<xsl:variable name="ctype" select="@complex_type"/>
	<xsl:choose>
		<xsl:when test="@ecl_type"><xsl:value-of select="@ecl_type"/></xsl:when>
		<xsl:otherwise>
			<xsl:variable name="srcfile" select="/expesdl/types/type[@name=$ctype]/@src"/>
			<xsl:if test="$sourceFileName != $srcfile">
				<xsl:value-of select="$srcfile"/><xsl:text disable-output-escaping="yes">.</xsl:text>
			</xsl:if>
			<xsl:text disable-output-escaping="yes">t_</xsl:text><xsl:value-of select="$ctype"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="output_ecl_array_type">
	<xsl:variable name="ctype" select="@type"/>
	<xsl:choose>
		<xsl:when test="@ecl_type"><xsl:value-of select="@ecl_type"/></xsl:when>
		<xsl:otherwise>
			<xsl:variable name="srcfile" select="/expesdl/types/type[@name=$ctype]/@src"/>
			<xsl:if test="$sourceFileName != $srcfile">
				<xsl:value-of select="$srcfile"/><xsl:text disable-output-escaping="yes">.</xsl:text>
			</xsl:if>
			<xsl:text disable-output-escaping="yes">t_</xsl:text><xsl:value-of select="$ctype"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="output_comments">
	<xsl:if test="@ecl_comment"><xsl:value-of select="@ecl_comment"/></xsl:if>
	<xsl:choose>
		<xsl:when test="@complex_type">
			<xsl:variable name="ctype" select="@complex_type"/>
			<xsl:if test="/expesdl/types/type[@name=$ctype]/@comment">
				<xsl:text disable-output-escaping="yes">//</xsl:text><xsl:value-of select="/expesdl/types/type[@name=$ctype]/@comment"/>
			</xsl:if>
		</xsl:when>
		<xsl:when test="@enum_type">
			<xsl:variable name="etype" select="@enum_type"/>
			<xsl:if test="/expesdl/types/type[@name=$etype]/@comment">
				<xsl:text disable-output-escaping="yes">//</xsl:text><xsl:value-of select="/expesdl/types/type[@name=$etype]/@comment"/>
			</xsl:if>
		</xsl:when>
	</xsl:choose>
	<xsl:if test="@optional">
		<xsl:text disable-output-escaping="yes">//hidden[</xsl:text><xsl:value-of select="@optional"/><xsl:text disable-output-escaping="yes">]</xsl:text>
	</xsl:if>
	<xsl:if test="@ecl_type and (@type or @complex_type)">
        <xsl:choose>
         <xsl:when test="name()='EsdlArray|EsdlList'">
           <xsl:text disable-output-escaping="yes"> // Real type: </xsl:text>
           <xsl:text disable-output-escaping="yes">DATASET(t_</xsl:text>
           <xsl:value-of select="@type"/>
           <xsl:text disable-output-escaping="yes">) </xsl:text>
           <xsl:call-template name="output_ecl_name"/>
           <xsl:text disable-output-escaping="yes"> {XPATH('</xsl:text>
           <xsl:choose>
		<xsl:when test="name() ='EsdlArray'">
                  <xsl:if test="not(@flat_array)"><xsl:value-of select="@name"/></xsl:if>
                  <xsl:text disable-output-escaping="yes">/</xsl:text>
                  <xsl:call-template name="output_item_tag"/>
                </xsl:when>
		<xsl:otherwise>
                  <xsl:value-of select="@name"/>
	        </xsl:otherwise>
             </xsl:choose>
           <xsl:text disable-output-escaping="yes">')};</xsl:text>
         </xsl:when>
         <xsl:when test="name()='EsdlElement' and starts-with(@ecl_type,'tns:')">
           <xsl:text disable-output-escaping="yes"> // Real type: RECORD t_</xsl:text>
           <xsl:value-of select="@type|@complex_type"/>
         </xsl:when>
         <xsl:when test="name()='EsdlElement' and not(starts-with(@ecl_type,'tns:'))">
           <xsl:text disable-output-escaping="yes"> // Xsd type: </xsl:text>
           <xsl:value-of select="@type|@complex_type"/>
         </xsl:when>
         <xsl:otherwise>
             <xsl:value-of select="@type"/>
         </xsl:otherwise>
        </xsl:choose>
    </xsl:if>
</xsl:template>

<xsl:template name="output_xpath">
    <xsl:text disable-output-escaping="yes">XPATH('</xsl:text>
	<xsl:choose>
		<xsl:when test="@ecl_path"><xsl:value-of select="@ecl_path"/></xsl:when>
		<xsl:otherwise><xsl:if test="@attribute"><xsl:value-of select="'@'"/></xsl:if> <xsl:value-of select="@name"/></xsl:otherwise>
	</xsl:choose>
	<xsl:text disable-output-escaping="yes">')</xsl:text>
</xsl:template>

<xsl:template name="output_item_tag">
    <xsl:choose>
         <xsl:when test="@item_tag"><xsl:value-of select="@item_tag"/></xsl:when>
         <xsl:when test="@type and (@type='int' or @type='integer' or @type='bool' or @type='short' or @type='float' or @type='double' or @type='string' or @type='long' or @type='decimal' or @type='byte' or @type='unsignedInt' or @type='unsignedShort' or @type='unsignedByte')">
           <xsl:value-of select="'Item'"/>
         </xsl:when>
         <xsl:when test="@type">
           <xsl:value-of select="@type"/>
         </xsl:when>
         <xsl:otherwise><xsl:value-of select="'Item'"/></xsl:otherwise>
    </xsl:choose>
</xsl:template>

</xsl:stylesheet>
